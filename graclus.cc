#include <assert.h>
#include "tensorflow/core/framework/op.h"
#include "tensorflow/core/framework/op_kernel.h"
#include "tensorflow/core/framework/shape_inference.h"
#include "tensorflow/core/util/work_sharder.h"
#include "tensorflow/core/util/tensor_format.h"
#include <iostream>
#include <metis.h>

using std::vector;
using namespace tensorflow;
using shape_inference::Shape;
using shape_inference::Dimension;
using shape_inference::DimensionHandle;
using shape_inference::ShapeHandle;

int boundary_points;
int spectral_initialization;
int cutType; /*cut type, default is normalized cut */
int memory_saving; /* forbid using local search or empty cluster removing */
/*char mlwkkm_fname[256]; */ /*used to store coarsest file*/

/*************************************************************************
* multi-level weighted kernel k-means main function
**************************************************************************/
idxtype* graclus(std::vector<int>& xadj, std::vector<int>& adjncy, std::vector<int>& adjwgt, int n_clust)
{
  int edgenum, curredge, numedges, currval;

  int readew, readvw, fmt, ncon, tmp;
  int i, nparts=1, options[11];
  idxtype *part;
  float rubvec[MAXNCON], lbvec[MAXNCON];
  float result;
  GraphType graph;
  int numflag = 0, wgtflag = 0, edgecut, chain_length;
  int no_args = 1, clusteringEva =0, levels;
  idxtype *t1, *t2;

  edgenum = adjncy.size();
  nparts = n_clust;
  cutType = 0;
  chain_length = 0;
  spectral_initialization = 0;
  
  memory_saving = 0;
  boundary_points = 0;

  /*printf("edgenum = %d, edgenum2 = %d\n", edgenum, mxGetNumberOfElements(prhs[0]));*/
   readew=1;
   InitGraph(&graph);

   graph.nvtxs = adjncy.size();
   graph.nedges = edgenum;
   graph.xadj = xadj.data();
   graph.adjncy = adjncy.data();
   graph.adjwgt = adjwgt.data();

   ncon = 0;
   fmt = 1;
   ncon = graph.ncon = (ncon == 0 ? 1 : ncon);
   readew = (fmt%10 > 0);
   readvw = ((fmt/10)%10 > 0);
   wgtflag = 0;
   if (readew)
     wgtflag += 1;
   if (readvw)
     wgtflag += 2;

   /*levels = 5*nparts;*/
   levels = amax((graph.nvtxs)/(40*log2_metis(nparts)), 5*(nparts));
   /*printf("Will coarsen until %d nodes...\n", levels);*/

  part = (idxtype*) malloc(sizeof(idxtype)*(graph.nvtxs));
  options[0] = 0;
    if (graph.ncon == 1) 
    {
      /*METIS_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
	&wgtflag, &numflag, &nparts, options, &edgecut, part); 
      */
      MLKKM_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
			  &wgtflag, &numflag, &nparts, &chain_length, options, &edgecut, part, levels);
    
      /* ends */
    }
    else {
      for (i=0; i<graph.ncon; i++)
	rubvec[i] = HORIZONTAL_IMBALANCE;
      /*
	METIS_mCPartGraphKway(&graph.nvtxs, &graph.ncon, graph.xadj, graph.adjncy, graph.vwgt, 
	graph.adjwgt, &wgtflag, &numflag, &nparts, rubvec, options, &edgecut, part);
      */
    }
  ComputePartitionBalance(&graph, nparts, part, lbvec);
  if (cutType == 0){
    result = ComputeNCut(&graph, part, nparts);
  }
  else{
    //result = ComputeRAsso(&graph, part, nparts);
  }
  return part;
}  


REGISTER_OP("Graclus")
    .Input("adj: bool")
    .Input("weights: float32")
    .Input("n_clusters: int64")
    .Output("supernode_assign: int32")
    .SetShapeFn([](::tensorflow::shape_inference::InferenceContext* c) {        
        c->set_output(0, c->input(0));
        return Status::OK();
    });

class GraclusOp : public OpKernel {
 public:
      explicit GraclusOp(OpKernelConstruction* context) : OpKernel(context) { }

  void Compute(OpKernelContext* context) override {
    // Grab the input tensor
    const Tensor& adj = context->input(0);
    auto adj_arr = adj.tensor<bool, 3>();
    const Tensor& weights = context->input(1);
    auto weights_arr = weights.tensor<float, 3>();
    const Tensor& n_clusters = context->input(2);
    auto n_clusters_arr = n_clusters.tensor<int, 1>();
    // Create an output tensor
    Tensor* supernode_assign = NULL;
      
    const TensorShape& adj_shape = adj.shape();
    TensorShape out_shape = adj_shape;
    OP_REQUIRES_OK(context, context->allocate_output(0, out_shape,
                                                     &supernode_assign));
    
    auto output_shape = out_shape.dim_sizes();
    auto output_flat = supernode_assign->flat<int>();

    // Set all but the first element of the output tensor to 0.
    const int N = output_flat.size();
    for (int i = 1; i < N; i++) {
      output_flat(i) = 0;
    }

    auto supernode_assign_arr = supernode_assign->tensor<int, 3>();

    const int batch_size = output_shape[0];
    std::function<void(int64, int64)> shard;
    shard = [&adj_arr, &weights_arr, &n_clusters_arr, &supernode_assign_arr, &output_shape](int64 start, int64 limit) {
        for (int graph = start; graph < limit; ++graph) {
            int n_clus = n_clusters_arr(graph);
            std::vector<int> xadj;
            xadj.push_back(0);
            std::vector<int> adjncy;
            std::vector<int> adjwgt;
            for (int m = 0; m < output_shape[1]; m++) {
                for (int n = 0; n < output_shape[0]; n++) {
                    if (m != n && adj_arr(graph, m, n)) {
                        adjncy.push_back(n);
                        adjwgt.push_back(weights_arr(graph, m, n));
                    }
                }
                xadj.push_back(xadj.size());
            }
            idxtype* assign = graclus(xadj, adjncy, adjwgt, n_clus);
            
            for(int i = 0; i < adjncy.size(); i++) {
                supernode_assign_arr(graph, i, assign[i]) = 1.;
            }
            free(assign);
        }
    };

    // This is just a very crude approximation
    const int64 single_cost = 10000 * output_shape[1] * output_shape[2];

    auto worker_threads = context->device()->tensorflow_cpu_worker_threads();
    Shard(worker_threads->num_threads, worker_threads->workers, batch_size, single_cost, shard);
  }

 private:
    TensorFormat data_format_;
};

REGISTER_KERNEL_BUILDER(Name("Graclus").Device(DEVICE_CPU), GraclusOp);

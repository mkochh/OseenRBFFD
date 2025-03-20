#include  "../../H2Lib/h2lib.h"

#ifndef OSEENCLUSTER_HEADER

/**
 * \brief Inline utility function for setting an auxilary array for domain decomposition clustering
 * 
 * \param idx Index set to set the flags for.
 * \param size Size of the index set.
 * \param sp Sparse matrix defining connections. 
 * \param flag Auxilary array to be set.
 */
inline void set_flags(uint *idx, uint size, psparsematrix sp, uint *flag);

/**
 * \brief Reset an auxilary array for domain decomposition based clustering
 *  All entries are set to zero.
 * 
 * \param flag_size Size of the auxilary array.
 * \param flag Auxilary array.
 */
inline void reset_flag(uint flag_size, uint *flag);

/**
 * \brief Auxilary function for swapping two indices at given 
 *  adresses.
 * 
 * \param i First index.
 * \param j Second index.
 */
inline void swap_indices(uint *i, uint *j);

/**
 * \deprecated Deprecated cluster function. Use \ref build_adaptive_from_cluster instead.
 * \brief Build cluster trees for velocity and pressure discretizations
 *        in a coupled way.
 *
 * \param cgu Geometric information about the velocity discretization.
 * \param cgp Geometric information about the pressure discretization.
 * \param sizeu Size of the velocity cluster.
 * \param sizep Size of the pressure cluster.
 * \param idxu Velocity discretization index set.
 * \param idxp Pressure discretization index set.
 * \param clf Maximal size of leaf clusters.
 * \param sp System matrix for the domain decomposition based clustering
 *           of the veloctity discretizations index set.
 * \param dim Spatial dimension of the discretization.
 * \param flag Auxilary array for the domain decomposition based clustering.
 * \return pcluster* The velocity and the pressure cluster tree.
 *
 */
pcluster *build_coupled_adaptive_dd_cluster(pclustergeometry cgu,
                                            pclustergeometry cgp, uint sizeu,
                                            uint sizep, uint *idxu, uint *idxp,
                                            uint clf, psparsematrix sp,
                                            uint dim, uint *flag);

/**
 * \brief Build a \ref cluster tree with a domain decomposition based method.
 *        The interface clusters are build broader compared to the standard
 *        variant by choosing nodes for the interface from both domain clusters.
 *
 * \param cg \ref clustergeometry object containing geometric informations.
 * \param size Size of the cluster.
 * \param idx Index set of the cluster.
 * \param clf Maximal leaf size.
 * \param sp System matrix determining connection between the nodes of the cluster.
 * \param flag Auxillary array for building the interface cluster.
 * \return pcluster \ref cluster tree for the index set.
 */
pcluster build_adaptive_broad_dd_cluster(pclustergeometry cg, uint size,
                                         uint *idx, uint clf, psparsematrix sp,
                                         uint *flag);

/* ********************************************************
 * Simultaniously subdivision based clustering methods    *
 **********************************************************/

/**
 * \brief Build a \ref cluster tree with a domain decomposition based method.
 *        The sons are splitted simultaniously, so that the splitting height
 *        is the same for all sons.
 *
 * \param cg \ref clustergeometry object containing geometric informations.
 * \param size Size of the cluster.
 * \param idx Index set of the cluster.
 * \param clf Maximal leaf size for the cluster tree.
 * \param sp System matrix determining connection between the nodes of the cluster.
 * \param direction Splitting direction.
 * \param flag Auxilary array.
 * \return pcluster \ref cluster tree for the index set.
 */
pcluster build_simsub_dd_cluster(pclustergeometry cg, uint size,
                                 uint *idx, uint clf, psparsematrix sp,
                                 uint direction, uint *flag);

/**
 * \brief Build a \ref cluster tree with simultanious sub division. Compared to the
 *        H2Lib method, the level is increased with every division in one direction,
 *        instead after divisions in each direction.
 *
 * \param cg \ref clustergeometry object containing geometric informations.
 * \param size Size of the cluster.
 * \param idx Index set of the cluster.
 * \param clf Maximal leaf size.
 * \param direction Direction of partition.
 * \return pcluster \ref cluster tree for the index set.
 */
pcluster build_simsub_cluster_level(pclustergeometry cg, uint size,
                                    uint *idx, uint clf, uint direction);

pcluster build_coupled_dd_cluster(pcluster cluster_p, pclustergeometry cg, uint size,
                                  uint *idx, uint clf, psparsematrix graph, uint *flag);

/**
 * \brief Build a \ref cluster tree coupeld with a cluster tree build with bisection. 
 * 
 * This is the second version of the coupled clustering. The index set is subdivided 
 * into five sub clusters w.r.t. a second cluster tree for a different index set.
 * The coupling is based on connections given by a sparse matrix. The interface 
 * is split into three parts, based on geometric informations given by a clustergeometry 
 * object and another sparse matrix.  
 * 
 * \param cluster_p Associated (coupled) cluster tree. Has to be built by bisection previously.
 * \param cg Clustergeometry object containing geometric informations. 
 * \param size Size of the index set.
 * \param idx Index set.
 * \param clf Maximal leaf size for the interface clustering.
 * \param sp_coupling Sparse matrix indicating connections between the 
 * index set and the associated cluster.
 * \param sp_vel Sparse matrix indicating connections between the indices.
 * \param flag Auxilary array.
 * \return pcluster \ref cluster tree for the index set.
 */
pcluster build_coupled_dd_cluster_v2(pcluster cluster_p, pclustergeometry cg, uint size, uint *idx, 
                                     uint clf, psparsematrix sp_coupling, psparsematrix sp_vel, uint *flag);

#ifdef USE_NETCDF
void write_cdf_cluster_grid(pccluster t, pclustergeometry cg, const char *name);
#endif
#endif
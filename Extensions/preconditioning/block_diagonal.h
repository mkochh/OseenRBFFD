#include "prec.h"

/**
 * @brief Class representing block diagonal preconditioners for saddle point
 *        matrices, where the approximates for the upper left and the
 *        Schur complement block are given as preconditioners of type
 * Preconditioner
 *
 */
class Block_Diagonal_Prcd : public Preconditioner
{
public:
  Preconditioner *prcd_A;
  Preconditioner *prcd_S;

  uint n; // Number of rows/columns for the upper left block
  uint m; // Number of rows/columns for the Schurcomplement block

  /**
   * @brief Construct a new block diagonal preconditioner for saddle point
   * systems
   *
   * @param prec_A Preconditioner for the upper left block
   * @param prec_S Preconditioner for the Schurcomplement
   */
  Block_Diagonal_Prcd(Preconditioner *prcd_A, uint n, Preconditioner *prcd_S,
                      uint m);
  /**
   * @brief Destroy the Block_Diagonal_Prec object
   *
   */
  ~Block_Diagonal_Prcd();

  /**
   * @brief Apply the block preconditioner to a given vector
   *
   * @param r Vector to apply the preconditioner on
   */
  void apply_preconditioner(pavector r);
};
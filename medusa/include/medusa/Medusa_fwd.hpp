#ifndef MEDUSA_MEDUSA_FWD_HPP_
#define MEDUSA_MEDUSA_FWD_HPP_

/**
 * @file
 * Main header file that includes only declarations for some heavy classes such as @ref mm::WLS
 * to reduce compile time.
 *
 * Linker errors are possible when using only this header due to missing definitions of some
 * classes and functions. There are two solutions for this, either include full header
 * file Medusa.hpp or explicitly instantiate missing classes in appropriate `.cpp` files.
 */

#include "Config.hpp"
#include "Domain_fwd.hpp"
#include "Approximations_fwd.hpp"
#include "Interpolants_fwd.hpp"
#include "Types.hpp"
#include "Integrators.hpp"
#include "IO.hpp"
#include "Operators.hpp"
#include "Utils.hpp"

#endif  // MEDUSA_MEDUSA_FWD_HPP_

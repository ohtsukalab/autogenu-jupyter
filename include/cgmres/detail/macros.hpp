#ifndef CGMRES__MACROS_HPP_
#define CGMRES__MACROS_HPP_

namespace cgmres {
namespace detail {

#define CGMRES_EIGEN_CONST_CAST(TYPE, OBJ) const_cast<TYPE &>(OBJ.derived())

} // namespace detail
} // namespace cgmres

#endif // CGMRES__MACROS_HPP_
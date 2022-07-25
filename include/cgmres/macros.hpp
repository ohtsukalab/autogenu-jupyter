#ifndef MACROS_HPP_
#define MACROS_HPP_

namespace cgmres {

#define EIGEN_CONST_CAST(TYPE, OBJ) const_cast<TYPE &>(OBJ.derived())

} // namespace cgmres

#endif // MACROS_HPP_
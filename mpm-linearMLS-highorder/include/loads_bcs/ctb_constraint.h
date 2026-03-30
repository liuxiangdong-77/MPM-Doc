#ifndef MPM_CTB_CONSTRAINT_H_
#define MPM_CTB_CONSTRAINT_H_

#include "data_types.h"

namespace mpm {

//! CTBConstraint class to store friction constraint on a set
//! \brief CTBConstraint class to store a constraint on a set
//! \details CTBConstraint stores the constraint as a static value
//! 负责接收json输入文件参数的类
class CTBConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] dir Direction of p-wave propagation in model
  //! \param[in] n 透射次数
  //! \param[in] position Nodal position along boundary
  CTBConstraint(int setid, unsigned dir, unsigned int order, double delta, double h_min,
                mpm::Edge_Position position = mpm::Edge_Position::None,
                mpm::Position position_simple = mpm::Position::None)
      : setid_{setid},
        dir_{dir},
        order_{order},
        delta_{delta},
        h_min_{h_min},
        position_{position},
        position_simple_ {position_simple} {};

  // Set id
  int setid() const { return setid_; }

  // Direction of p-wave travel
  unsigned dir() const { return dir_; }

  //N
  unsigned int order() const { return order_; }

  // Virtual viscous layer thickness
  double delta() const { return delta_; }

  // Cell height
  double h_min() const { return h_min_; }

  // Return position
  mpm::Edge_Position position() const { return position_; }

  mpm::Position position_simple() const { return position_simple_;}

 private:
  // ID
  int setid_;
  // Direction
  unsigned dir_;
  // N
  unsigned int order_; //透射次数

  double delta_;

  double h_min_;

  // Node position
  mpm::Edge_Position position_;

  mpm::Position position_simple_;
};
}  // namespace mpm
#endif  // MPM_CTB_CONSTRAINT_H_ //lxd
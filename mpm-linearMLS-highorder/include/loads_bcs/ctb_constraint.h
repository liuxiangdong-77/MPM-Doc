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
  //! \param[in] n 透射次数
  //! \param[in] position Nodal position along boundary
  CTBConstraint(int setid, unsigned int order,
                mpm::Edge_Position position = mpm::Edge_Position::None)
      : setid_{setid},
        order_{order},
        position_{position} {};

  // Set id
  int setid() const { return setid_; }

  //N
  unsigned int order() const { return order_; }

  // Return position
  mpm::Edge_Position position() const { return position_; }

 private:
  // ID
  int setid_;
  // N
  unsigned int order_; //透射次数

  // Node position
  mpm::Edge_Position position_;
};
}  // namespace mpm
#endif  // MPM_CTB_CONSTRAINT_H_ //lxd
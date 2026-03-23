#ifndef MPM_MATERIAL_LINEAR_ELASTIC_CTB_H_
#define MPM_MATERIAL_LINEAR_ELASTIC_CTB_H_

#include <limits>
#include <vector>

#include "Eigen/Dense"

#include "linear_elastic.h"

namespace mpm {

//! LinearElasticCTB class for Characteristic Transmitting Boundary (CTB)
//! \brief Linear Elastic material model for Characteristic Transmitting Boundary
//! \details LinearElasticCTB class stores multiple time steps of state variables
//! for CTB implementation, inherits all mechanical behavior from LinearElastic
//! \tparam Tdim Dimension
template <unsigned Tdim>
class LinearElasticCTB : public LinearElastic<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;
  //! Define a vector of dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id
  //! \param[in] material_properties Material properties
  LinearElasticCTB(unsigned id, const Json& material_properties);

  //! Destructor
  ~LinearElasticCTB() override{};

  //! Delete copy constructor
  LinearElasticCTB(const LinearElasticCTB&) = delete;

  //! Delete assignement operator
  LinearElasticCTB& operator=(const LinearElasticCTB&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! State variables
  std::vector<std::string> state_variables() const override;
  
  unsigned history_size() const { return history_size_; }

  //! Advance history variables (call this at the end of each time step)
  //  预制历史变量（在每个时间步结束时调用）
  //! \param[in] state_vars Current state variables
  //! \param[in] particle Particle pointer to get current physical quantities
  void advance_history_variables(mpm::dense_map* state_vars, 
                                const ParticleBase<Tdim>* particle,
                                unsigned long long time_step);

  //! Get history variable at specific time step
  //  获取特定时间步的历史变量
  //! \param[in] state_vars State variables
  //! \param[in] var_name Variable name
  //! \param[in] time_step Time step offset (0 = current, 1 = previous, etc.)
  //! \retval value Variable value at specified time step
  double get_history_variable(const mpm::dense_map* state_vars,
                             const std::string& var_name, 
                             unsigned time_step) const;

  //! Check if material supports CTB history
  bool supports_ctb_history() const override { return true; }

 protected:
  //! material id
  using LinearElastic<Tdim>::id_;
  //! Material properties
  using LinearElastic<Tdim>::properties_;
  //! Logger
  using LinearElastic<Tdim>::console_;

 private:
  //! Number of history steps to store
  unsigned int history_size_{};//constexpr: 编译时确定的常量表达式
};  // LinearElasticCTB class
}  // namespace mpm

#include "linear_elastic_ctb.tcc"

#endif  // MPM_MATERIAL_LINEAR_ELASTIC_CTB_H_
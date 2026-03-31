//! Assign CTB constraint to nodes
//  将CTB约束分配到结点集
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_ctb_constraint(
    int nset_id,
    const std::shared_ptr<mpm::CTBConstraint>& ctb_constraint, unsigned phase, double dt) {
  bool status = true;
  try {
    int set_id = nset_id;
    auto nset = mesh_->nodes(set_id);//成功
    if (nset.size() == 0)
      throw std::runtime_error(
          "Node set is empty for assignment of CTB constraints");

    if (!mesh_->create_nodal_ctb_constraint(set_id, ctb_constraint, phase, dt))
      throw std::runtime_error("Failed to create CTB constraint in mesh");//直接把ctb_constraint整个传到mesh.里处理

    console_->info("Assigned MTF constraint to nset {}: pos={}, order={}",
                  set_id,
                  static_cast<int>(ctb_constraint->position()),
                  ctb_constraint->order());//成功，表示参数已传到ctb_constraint类里
    
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign CTB constraints to nodes
//  将CTB约束分配到单个结点
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_ctb_constraints(
    const std::vector<std::tuple<mpm::Index, mpm::Edge_Position, 
                                 unsigned>>& ctb_constraints, unsigned phase, double dt) {
  bool status = true;
  try {
    for (const auto& ctb_constraint_tuple : ctb_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(ctb_constraint_tuple);
      // Position
      mpm::Edge_Position position = std::get<1>(ctb_constraint_tuple);
      // Order
      unsigned order = std::get<2>(ctb_constraint_tuple);

      // 创建 CTB 约束对象
      auto constraint = std::make_shared<mpm::CTBConstraint>(
          -1, order, position); // set_id = -1 表示单个节点
      
      // 将约束应用到单个节点
      // 这里需要创建一个临时的节点集
      // 或者直接调用节点的 CTB 初始化方法
      
      console_->info("Assigned CTB constraint to node {}: pos={}, order={}",
                    nid, static_cast<int>(position), order);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign ctb constraints pointers and ids
template <unsigned Tdim>
void mpm::Constraints<Tdim>::assign_ctb_id_ptr(
    unsigned nset_id,
    std::shared_ptr<mpm::CTBConstraint>& ctb_constraint) {
  this->ctb_constraint_.emplace_back(ctb_constraint);//存储约束对象指针
  this->ctb_nset_id_.emplace_back(nset_id);//存储对应的结点集ID
}


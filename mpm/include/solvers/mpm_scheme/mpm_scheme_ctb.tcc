//! Store particle history for CTB (根据order决定存储步数)
template <unsigned Tdim>
inline void mpm::MPMScheme<Tdim>::store_particle_history_for_ctb(unsigned long long time_step) {
  
  // 遍历所有粒子，存储历史数据
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::advance_ctb_history,
                std::placeholders::_1, time_step));

}
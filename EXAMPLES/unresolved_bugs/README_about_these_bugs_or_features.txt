
- bugPMLs_1 : Strong spurious reflections from the PMLs in some particular cases.

- bugPMLs_2 : Weak instability slowly developing in long times. This problem (the worst by far) occurs in simulations with many time steps (> 100000) and viscoelastic media with PMLs. Even if the PMLs are not in contact with the viscoelastic part!

- bugViscoelasticity : This shows a gradual, very slow explosion in viscoelastic media. This issue is probably not a "bug" strictly speaking but may arise from the inconsistency of memory variables between two adjacent elements. The real problem occurs when considering also acoustic media in the same run. In that case the instability occurs faster.


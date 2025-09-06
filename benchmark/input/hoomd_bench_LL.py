# hoomd_benchmark_cpu1.py
# LJ benchmark matching your C++: 20x20x20 SC lattice, mass=1, v=0,
# LJ(eps=5, sig=1) truncated (unshifted) at r_cut=3*sig,
# Cell neighbor list buffer=0, ConstantVolume (NVE), dt=2e-4, steps=10000.
# CPU-only, single thread.

import os
# Force one thread (HOOMD respects OMP_NUM_THREADS on CPU)
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import time
import numpy as np
import hoomd
import hoomd.md


for n in [[20,20,20], [20,15,15], [10,15,15], [10,10,10]]:
    print()
    for _ in range(5):
        # ----------------------------
        # Parameters (match your C++)
        # ----------------------------
        NX, NY, NZ = n
        a = 1.1225
        sigma = 1.0
        epsilon = 5.0
        r_cut = 3.0 * sigma

        dt = 0.0002
        steps = 10_000

        # Grid physical span: (N-1)*a
        Lx = (NX - 1) * a
        Ly = (NY - 1) * a
        Lz = (NZ - 1) * a

        # Box with margin >= r_cut around grid (your code uses 1.5x)
        extent_scale = 2
        Ex = extent_scale * Lx
        Ey = extent_scale * Ly
        Ez = extent_scale * Lz

        # ----------------------------
        # Device (CPU)
        # ----------------------------
        device = hoomd.device.CPU(notice_level=1)
        sim = hoomd.Simulation(device=device, seed=42)

        # ----------------------------
        # Build simple-cubic lattice
        # ----------------------------
        positions = np.empty((NX * NY * NZ, 3), dtype=np.float32)
        idx = 0
        x0, y0, z0 = -0.5 * Lx, -0.5 * Ly, -0.5 * Lz
        for i in range(NX):
            for j in range(NY):
                for k in range(NZ):
                    positions[idx] = (x0 + i * a, y0 + j * a, z0 + k * a)
                    idx += 1
        N = positions.shape[0]

        # Initialize from a Snapshot (set box + particle data on rank 0)
        snap = hoomd.Snapshot()
        if snap.communicator.rank == 0:
            snap.configuration.box = [Ex, Ey, Ez, 0.0, 0.0, 0.0]  # [Lx, Ly, Lz, xy, xz, yz]
            snap.particles.N = N
            snap.particles.types = ['A']
            snap.particles.typeid[:] = np.zeros(N, dtype=np.uint32)
            snap.particles.position[:] = positions
            snap.particles.mass[:] = np.ones(N, dtype=np.float32)
            snap.particles.velocity[:] = np.zeros((N, 3), dtype=np.float32)

        sim.create_state_from_snapshot(snapshot=snap)

        # ----------------------------
        # Forces: Lennard-Jones 12-6
        # ----------------------------
        # Cell neighbor list with buffer=0.0 (rebuilds as needed ~every step for tiny dt)
        nl = hoomd.md.nlist.Cell(buffer=0.3)  # for more performance skin = 0.3-0.5
        lj = hoomd.md.pair.LJ(nlist=nl, mode='none')  # 'none' = unshifted, plain truncation
        lj.params[('A', 'A')] = dict(epsilon=epsilon, sigma=sigma)
        lj.r_cut[('A', 'A')] = r_cut

        # ----------------------------
        # Integrator: NVE (ConstantVolume)
        # ----------------------------
        integrator = hoomd.md.Integrator(dt=dt)
        method = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())  # NVE (no thermostat)
        integrator.methods = [method]
        integrator.forces = [lj]
        sim.operations.integrator = integrator

        # ----------------------------
        # Run benchmark & report
        # ----------------------------
        t0 = time.perf_counter()
        sim.run(steps)
        elapsed = time.perf_counter() - t0
        tps = steps / elapsed if elapsed > 0 else float('inf')

        print(f"Ran {steps} steps for N={N} in {elapsed:.3f} s => {tps:,.0f} steps/s")
        # print(f"Device: CPU | OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}")
        print(f"Neighbor list buffer: {nl.buffer} | LJ mode: {lj.mode} | r_cut={r_cut}")

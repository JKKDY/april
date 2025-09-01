# APRIL -  A Particle Runtime Interaction Library

APRIL is a small, modular C++ library for particle-based simulations.
It aims to combine high performance with a flexible, easy-to-use, and expressive API. The library emphasizes clear architecture, plug-and-play components, and modern C++ features (concepts, CRTP-style dispatch).
> Status: Most of the foundation in place (environment, interactions, containers, integrator, monitors). Additional foundational features are in development (controllers, force fields, boundary conditions)


## Core Features

- **Modular design**: swap or extend **forces**, **containers** (force calculators), **integrators**, and **monitors**.
- **Modern C++**: concepts for compile-time interface checking; a simple, readable public API via `april/april.h`.
- **Ergonomic setup**: clear setup path. Special care was taken to minimize template verbosity with CTAD.
- **Built-in monitors**: binary snapshots, terminal diagnostics, progress bar, and a simple benchmark.
- **Built-in containers**: `DirectSum` (all-pairs) and `LinkedCells` (cell lists).
- **Tested core**: GoogleTest suite covering interactions, containers, integrator steps, binary I/O, and utilities.
- **Small animation script**: a Python helper to quickly preview simulation output.


## Getting Started

**Requirements:**
- C++20-capable compiler (e.g., GCC 12+, Clang 14+, MSVC 19.36+)
- CMake ≥ 3.20
- (Dev) GoogleTest for running the test suite

**Build Guide:**
From the project root run:
```bash
cmake -S . -B build
cmake --build build -j
```

To run tests use
```bash
ctest --test-dir build --output-on-failure
```

## Architecture

The following diagram shows the typical flow of a program using APRIL:
```
[Particles]   [Forces]          
     \          /                          
      v        v                           
   +-------------+            +-----------+ 
   | Environment |            | Container | 
   +-------------+            +-----------+ 
           \                       /         
            \                     /          
             v                   v           
            +---------------------+     transforms user-provided data (particles, 
            |   build_system(...) |  <- types/IDs, forces, domain) into dense internal
            +---------------------+     representations and wires the components together.
                       |
                       v
                  +---------+
                  | System  |   <— uses Container to compute forces
                  +---------+
                       ^
                       |  (each step: system.update_forces)
                       |
                +-------------+
                | Integrator  |
                +-------------+
                       |
                       |  emits records every N steps
                       v
                 +-----------+
                 | Monitors  |
                 +-----------+
```

- **Environment** \
Defines the simulation setup: particles, simulation domain (origin/extent), and the list of interactions (by **type pair** or **id pair**).
- **System** \
Materializes the environment: maps user IDs/types to dense internals, builds interaction tables, finalizes the domain, and wires everything together.
- **Container** \ 
Owns internal particle storage/indices and computes pairwise forces (e.g., **DirectSum**, **LinkedCells**). It’s injected with the interaction manager built by the system.
- **Integrator** \
Advances the state in time (e.g., Stoermer–Verlet). On each step it updates positions/velocities and asks the system to refresh forces.
- **Monitors** \
Optional observers invoked every step (or every *N* steps). They’re the primary way to emit output (binary frames, progress, timing, etc.).

## Usage (minimal example)
```c++
#include <april/april.h>
using namespace april;

// Simulation of a simple sun-planet system
int main() {
    // 1) Define an environment: particles + forces
    Environment env (forces<InverseSquare>);
    env.add({0,0,0}, {0,0,0}, 1.0, /*type*/0);         // Sun
    env.add({1,0,0}, {0,1,0}, 1e-3, /*type*/0);        // Planet
    env.add_force(InverseSquare(), to_type(0));        // gravity for type 0

    // 2) Choose a container (force calculator) and build a system
    auto system = build_system(env, DirectSum());

    // 3) Integrate with Stoermer–Verlet and attach monitors
    StoermerVerlet integrator(system, monitors<ProgressBar, Benchmark>);
    integrator.add_monitor(ProgressBar(50));           // updated every 50 steps            
    integrator.add_monitor(Benchmark());               // simple timing
    integrator.run_for(0.01, 10.0);                    // dt=0.01, T=10
}
```
Further examples can be found in `examples/`:
- Halley’s Comet - small N-body system with gravitational InverseSquare forces.
- Two-Body Collision -  MD-style collision with LennardJones interactions and explicit domain extents.

## Design Notes
- The entire public API is collected in april/april.h, so users normally only need a single include.
- Clear separation of user-facing vs. internal API: users work with declarative structs (e.g., `Environment`, `Particle`, container config structs) which are then consumed to build the internal representations (`impl::Particle`, container internals). This keeps the public API ergonomic and stable while allowing optimized internal implementations.
- Concepts enforce component interfaces at compile time (IsForce, IsMonitor, IsSystem), making extension points explicit.
- Environment → System → Integrator is the central workflow. Systems always delegate force evaluation to a chosen Container.
- Interactions can be specified by type pair or by particle id pair; missing cross-type entries are derived via a mix function.
- BinaryOutput writes a compact, versioned binary format (positions as float, plus type/id/state) suitable for lightweight analysis or visualization.
- Defaults (e.g., automatic domain extent/origin) are provided, but can always be overridden explicitly.


## Extending APRIL

APRIL’s components are designed to be easy to implement and drop in. In the following all nested namespaces inside april - aside from ::impl:: - are omitted for clarity. 

### Custom force

Implement the call operator and a `mix` rule (used to derive cross-type interactions) and provide a `cutoff_radius`:

```c++
struct MyForce {
    double cutoff_radius = -1.0; // negative = no cutoff
    vec3 operator()(const impl::Particle& a,
                    const impl::Particle& b,
                    const vec3& r_ab) const noexcept;
    MyForce mix(const MyForce& other) const noexcept;
};
```

Register in an environment:

```c++
Environment env (forces<MyForce>);
env.add_force(MyForce{...}, to_type(...));
```

### Custom container

Inherit from `impl::Container<Config, Env>`and provide 

- `build(const std::vector<Particle>&)`
- `calculate_forces()`
- accessors for particle/id/index ranges

Additionally provide a struct pointing to the containers type. This is passed into the `Config` template parameter as well as into  `build_system(...).`

```c++
template <class Env> class MyContainerImpl; // forward declaration

struct MyContainer { // User facing declaration
	template<class  Env> using impl = MyContainerImpl<Env>;
	... // optional user config data here
};

template <class Env>
class MyContainerImpl final : public Container<MyContainer, Env> {
    using Base = Container<MyContainerCfg, Env>;
    using Particle = typename Base::Particle;
public:
    using Base::Base; // forward config. To access config data use this->config.
	
    void build(const std::vector<Particle> & particles) { ... }
    void calculate_forces() { ... }
    
    // provide stable id-based access and sequential index-based access
    Particle& get_particle_by_id(ParticleID id) noexcept { /* ... */ }
    Particle& get_particle_by_index(size_t idx) noexcept { /* ... */ }

    size_t particle_count() const { ... }
    
    ...
};
```

Usage: 

```
auto system = build_system(env, MyContainer());
```

You can derive from `april::cont::impl::ContiguousContainer<Config, Env>` to reuse storage and id/index utilities. ContiguousContainer stores particles in a single contiguous vector.

### Custom integrator

Inherit from `core::impl::Integrator<System, MonitorPack<...>>` and provide `integration_step()`.

````c++
template<core::IsSystem Sys, class Pack>
class MyIntegrator;

template<core::IsSystem Sys, class... Ms>
class MyIntegrator<Sys, io::MonitorPack<Ms...>>
  : public impl::Integrator<Sys, MonitorPack<Ms...>> {
    using Base = impl::Integrator<Sys, MonitorPack<Ms...>>;
    using Base::sys; using Base::dt;

public:
    void integration_step() { ... }
};
````

### Custom monitor

Inherit from `io::Monitor` and implement `record(...)`. Optionally implement `before_step(...)` and `finalize()`.

```c++
class MyMonitor : public io::Monitor {
public:
    explicit MyMonitor(std::size_t every = 1) : Monitor(every) {}

    void record(std::size_t step, double time,
                const std::vector<impl::ParticleView>& particles) {
        // emit logs, write files, aggregate stats, etc.
    }

    // optional
    // void before_step(...);
    // void finalize();
};
```

Usage: 
````c++
StoermerVerlet integrator(system, io::monitors<MyMonitor>);
integrator.add_monitor(MyMonitor{10});  // call every 10 steps
````


## Roadmap

Planned additions (subject to change):
- [ ] Boundaries & boundary conditions 
- [ ] More integrators: 
  - [ ] Yoshida4
  - [ ] Boris Pusher
- [ ] Barnes–Hut container
- [ ] Controllers: particle sources/sinks, thermostats
- [ ] Force fields, including time-dependent fields
- [ ] Interacting geometries (e.g., box obstacles)
- [ ] SOA support
- [ ] Extendable particles via template parameter (e.g. add charge property)
- [ ] Parallelism
- [ ] C++ Modules


#pragma once


namespace april::core::impl {
	using namespace env;
    using namespace io;

	struct IntegratorCreateInfo {
		Environment * environment;
        
        OutputWriterBase * writer;
        OutputWriterBase * checkpoint_writer;
		Thermostat * thermostat;
		std::vector<env::ConstantForce> external_forces;
		std::vector<Observables*> observables;
	};

	class IntegratorBase {
	public:
		

		IntegratorBase(const IntegratorCreateInfo& create_info);

        /**
         * @brief Constructs a IntegratorBase object with a reference to a ParticleContainer and OutputWriter.
         * @param environment physical system to be simulated.
         * @param writer Pointer to an OutputWriterBase for writing output files.
         * @param checkpoint_writer Pointer to an CheckpointWriter to save a simulation state.
         * @param thermostat Thermostat to be used.
         * @param external_forces A list of external forces applied to the particles during the simulation.
         * @param stats Pointer to a Statistics object for statistics.
         */
        IntegratorBase(env::Environment& environment,
            std::unique_ptr<io::OutputWriterBase> writer = nullptr,
            std::unique_ptr<io::CheckpointWriter> checkpoint_writer = nullptr,
            const env::Thermostat& thermostat = env::Thermostat(),
            const std::vector<env::ConstantForce>& external_forces = {},
            std::unique_ptr<core::Statistics> stats = nullptr);
		
        virtual ~IntegratorBase() = default;
    protected:


	};
} // namespace april::core
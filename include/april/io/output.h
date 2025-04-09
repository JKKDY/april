#pragma once


namespace april::io {

	/**
	* @brief Abstract base class for output writers.
	* is used to generate output data during integration.
	*/
	class OutputWriterBase {
	public:
		
		/**
		* @brief Constructor for the OutputWriterBase class.
		* @param write_freq The frequency at which the output writer writes data.
		*/
		OutputWriterBase(const unsigned int write_freq) : write_frequency(write_freq) {}
		virtual ~OutputWriterBase() = default;

		/**
		* @brief creates output data. is called every <write_frequency> iterations
		* @param environment The current environment of the simulation.
		* @param iteration The current simulation iteration.
		*/
		virtual void write(const env::Environment& environment, int iteration) = 0;

		const unsigned int write_frequency; ///< The frequency at which the output writer writes data.
	};


	class CheckpointWriterBase {
	public:
		/**
		 * @brief Constructs a new CheckpointWriterBase object.
		 * @param file_base_name base name of the checkpoint files (default: "checkpoint").
		 * @param output_dir name of the output directory (default: "Checkpoint").
		 */
		explicit CheckpointWriterBase(std::string file_base_name = "checkpoint", std::string output_dir = "Checkpoints") :
			file_base_name(std::move(file_base_name)), output_dir(std::move(output_dir)) {}

		virtual ~CheckpointWriterBase() = default;

		/**
		 * @brief writes particle information to a checkpoint file.
		 * @param environment The current environment of the simulation.
		 * @param iteration The current simulation iteration.
		 */
		void write(const env::Environment& environment, int iteration) = 0;

	private:
		const std::string file_base_name;  ///< base name of the checkpoint files.
		const std::string output_dir;      ///< name of the output directory.
	};
} // namespace april::core
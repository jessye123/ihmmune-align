package iHMMuneAlign;

/**
 * @author harris
 *
 * TODO
 * 
 */
/**
 * container of static methods for finding the exponential decay in probability
 * for a nucleotide in a D or J gene
 */
public class ExponentialDecay {
	// main method for testing only
	public static void main(String[] args) {
		double VGene_nucl_prob = exponentialDecayVGene("acgacgacgacg", 10, 10);
		System.out.println("Dgene nucl prob = " + VGene_nucl_prob);
		
		double DGene_nucl_prob = exponentialDecayDGene("acgacgacgacg", 10, 100);
		System.out.println("Dgene nucl prob = " + DGene_nucl_prob);

		double JGene_nucl_prob = exponentialDecayJGene("acgacgacgacg", 10, 100);
		System.out.println("Jgene nucl prob = " + JGene_nucl_prob);

	}//--main()

	// the exponential decay rate found through data analyzis
	static final double EXPONENTIAL_DECAY_RATE = -0.0024;

	// since we do not know which DGene (it's length) nor the length
	// of the N region, we have to use an average expected length
	// for approximating the sequence position of DGene and JGene
	// nucleotides
	//
	// NB! for the VD-N region and the DJ-N region, we must take into
	// account the possibility of exonuclease removals at genes in both ends
	// For simplicity, I will assume that the Exo nucl. removal and the
	// N-addition cancel each other out
	// So we only have to find the average DGene length
	static final int average_DGene_length = 24;

	//	static final int average_Exonucleased_DGene_length = 15;
	//	static final int average_VD_N_Region_length = 4;
	//	static final int average_DJ_N_Region_length = 5;

	/**
	 * get the exponential decay in probability (e^(-0.0024*k) for a nucleotide
	 * in V-Gene "VGene" at position V-Gene position "VGene_nucl_position" with
	 * VGene start offset at "VGene_start_offset"
	 */
	public static double exponentialDecayVGene(String VGene,
			int VGene_nucl_position, int VGene_start_offset) {
		// find the sequence position absolute VGene nucl. position
		int sequence_position = VGene_nucl_position + VGene_start_offset;

		// get the exponential decay probability based on this position
		double exponential_decay_VGene_prob = exponentialDecay(sequence_position);

		// return this result
		return exponential_decay_VGene_prob;
	}//--exponential_decay_VGene_prob

	/**
	 * get the exponential decay in probability (e^(-0.0024*k) for a nucleotide
	 * in D-Gene "DGene" at position D-Gene position "DGene_nucl_position"
	 * knowing only the length of the VGene "k" can only be approximated because
	 * of our limited knowledge about the length of the N_Region(s)
	 */
	public static double exponentialDecayDGene(String DGene,
			int DGene_nucl_position, int completeVgeneLength) {
		// find the sequence position
		//int sequence_position = completeVgeneLength +
		// average_VD_N_Region_length + DGene_nucl_position;
		int sequence_position = completeVgeneLength + DGene_nucl_position;

		// get the exponential decay probability based on this position
		double exponential_decay_DGene_prob = exponentialDecay(sequence_position);

		// return this result
		return exponential_decay_DGene_prob;
	}//--exponential_decay_DGene_prob

	/**
	 * get the exponential decay in probability (e^(-0.0024*k) for a nucleotide
	 * in J-Gene "JGene" at position J-Gene position "DGene_nucl_position"
	 * knowing only the length of the VGene "k" can only be approximated because
	 * of our limited knowledge about the length of the N_Region(s)
	 */
	public static double exponentialDecayJGene(String JGene,
			int JGene_nucl_position, int completeVgeneLength) {
		// find the sequence position
		//	int sequence_position = completeVgeneLength +
		// average_VD_N_Region_length +
		//		average_DGene_length + average_DJ_N_Region_length +
		// JGene_nucl_position;
		int sequence_position = completeVgeneLength + average_DGene_length
				+ JGene_nucl_position;

		// get the exponential decay probability based on this position
		double exponential_decay_JGene_prob = exponentialDecay(sequence_position);

		// return this result
		return exponential_decay_JGene_prob;
	}//--exponential_decay_JGene_prob

	/**
	 * calculate e^(-0.0024*k) using "position_k" for k
	 */
	private static double exponentialDecay(int position_k) {
		double result = Math.exp(EXPONENTIAL_DECAY_RATE * position_k);

		// return the result found
		return result;

	}//--exponentialDecay()
}
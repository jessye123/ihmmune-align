package iHMMuneAlign;

/**
 * @author harris
 *
 * TODO
 * 
 */
/**
 * this interface is implemented by all classes delivering a mutability score
 */
interface MutabilityScore {
	/**
	 * get the mutability score for a nucleotide centered in the penta
	 * nucleotide "pentaNucleotide" returns the mutability score for the
	 * centered nucleotide
	 */
	public double pentaNucleotideScore(String pentaNucleotide);
}//--MutabilityScore

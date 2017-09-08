package iHMMuneAlign;

/**
 * @author harris
 *
 * HelperFunctions for finding codons
 * 
 */

/**
 * this class holds methods for calculating the mutability score of a
 * nucleotide, based on on whether it is part of a mutational hot spot or not.
 */
public class CodonHelperFunctions
{
	// /////////////////// Defines

	// //////// Nucleotides A,C,G,T
	final static char A_NUCLEOTIDE = 'a';

	final static char C_NUCLEOTIDE = 'c';

	final static char G_NUCLEOTIDE = 'g';

	final static char T_NUCLEOTIDE = 't';

	final static char N_NUCLEOTIDE = 'n';  // a/c/g/t
	
	final static int CODON_LENGTH = 3; // a codon is three nucleotides

	// ***********  amino acids
	public static boolean isAminoAcidW(char[] codon)
	{
		if(!isT(codon[0]))
			return false;
		if(!isG(codon[1]))
			return false;
		if(!isG(codon[2]))
			return false;
		
		return true;
	}//isAminoAcidW()
	
	public static boolean isAminoAcidG(char[] codon)
	{
		if(!isG(codon[0]))
			return false;
		if(!isG(codon[1]))
			return false;
		if(!isX(codon[2])) // last nucleotide can be A/C/G or T (X in other words)
			return false;
		
		return true;
	}//isAminoAcidW()
	
	// X is just any codon
	public static boolean isAminoAcidX(char[] codon)
	{
		if(!isX(codon[0]))
			return false;
		if(!isX(codon[1]))
			return false;
		if(!isX(codon[2]))
			return false;
		
		return true;
	}//isAminoAcidW()
	
	// ****** Stop Codon Functions
	
	public static boolean isTAG(char[] codon)
	{
		if(!isT(codon[0]))
			return false;
		if(!isA(codon[1]))
			return false;
		if(!isG(codon[2]))
			return false;
		
		return true;
	}//isTAG()
	
	public static boolean isTAA(char[] codon)
	{
		if(!isT(codon[0]))
			return false;
		if(!isA(codon[1]))
			return false;
		if(!isA(codon[2]))
			return false;
		
		return true;
	}//isTAA()
	
	public static boolean isTGA(char[] codon)
	{
		if(!isT(codon[0]))
			return false;
		if(!isG(codon[1]))
			return false;
		if(!isA(codon[2]))
			return false;
		
		return true;
	}//isTGA()
	
	
	// ****** single nucleotide functions
	
	/**
	 * test if nucleotide "nucleotide" is of type N, N_NUCLEOTIDE is only used
	 * at the end of a gene when a whole penta nucleotide is not available
	 * returns true or false
	 */
	private static boolean isN(char nucleotide) {
		if (nucleotide == N_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isA
	
	/**
	 * test if nucleotide "nucleotide" is of type A/C/G or T 
	 * returns true or false
	 */
	private static boolean isX(char nucleotide) {
		if (isA(nucleotide) || isC(nucleotide) || isG(nucleotide) || isT(nucleotide))
			return true;
		else return false;
	}//--isX

	/**
	 * test if nucleotide "nucleotide" is of type A returns true or false
	 */
	private static boolean isA(char nucleotide) {
		if(Character.toLowerCase(nucleotide) == A_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isA

	/**
	 * test if nucleotide "nucleotide" is of type C returns true or false
	 */
	private static boolean isC(char nucleotide) {
		if (Character.toLowerCase(nucleotide) == C_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isC

	/**
	 * test if nucleotide "nucleotide" is of type G returns true or false
	 */
	private static boolean isG(char nucleotide) {
		if (Character.toLowerCase(nucleotide) == G_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isG

	/**
	 * test if nucleotide "nucleotide" is of type T returns true or false
	 */
	private static boolean isT(char nucleotide) {
		if (Character.toLowerCase(nucleotide) == T_NUCLEOTIDE)
			return true;
		else
			return false;
	}//--isT

	// //////////////////////////////////////////////////////

	/**
	 * test if nucleotide "nucleotide" is of type R (A/G) returns true or false
	 */
	private static boolean isR(char nucleotide) {
		if (isA(nucleotide) || isG(nucleotide))
			return true;
		else
			return false;
	}//--isR

	/**
	 * test if nucleotide "nucleotide" is of type Y (C/T) returns true or false
	 */
	private static boolean isY(char nucleotide) {
		if (isC(nucleotide) || isT(nucleotide))
			return true;
		else
			return false;
	}//--isY

	/**
	 * test if nucleotide "nucleotide" is of type W (A/T) returns true or false
	 */
	private static boolean isW(char nucleotide) {
		if (isA(nucleotide) || isT(nucleotide))
			return true;
		else
			return false;
	}//--isW

}//--CodonHelperFunctions

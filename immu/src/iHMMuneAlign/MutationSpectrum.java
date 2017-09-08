package iHMMuneAlign;

/**
 * @author harris
 * 
 * TODO
 *  
 */
import java.io.*;
import java.util.*;

public class MutationSpectrum {

	public static void main(String[] args) {
		// create the mutationSpectrum class
		MutationSpectrum ms = new MutationSpectrum("C:\\Users\\test\\webworkspace\\immu\\src\\Mutation spectrum.txt");

		//	ms.getTNProbability(trinucleotide, mutateInto);

		/*
		 * // get the trinucleotide mutating String trinucleotide = args[0];
		 *  // get the symbol this Trinucleotide mutates into char mutateInto =
		 * args[1].charAt(0);
		 */

		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));

		String trinucleotide = null;
		String mutation = null;

		// keep reading from the keyboard until signaled to finish
		boolean finishedReading = false;
		while (!finishedReading) {
			System.out.println("enter the trinucleotide (\"END\" TO EXIT):  ");
			trinucleotide = getInputLine(br);

			if (trinucleotide.equalsIgnoreCase("end")) {
				finishedReading = true;
			} else {
				System.out.println("enter the mutation:  ");
				mutation = getInputLine(br);

				System.out.println();
				System.out.println();
				System.out.println("Searching for match");

				double prob = ms.getTNProbability(trinucleotide, mutation
						.charAt(0));

				System.out.println();
				System.out.println("finished Searching for match");
				System.out.println("prob = " + prob);
				System.out.println();
			}//--else

		}//--while
	}//--main() for testing only

	/**
	 * read a line from the keyboard
	 */
	private static String getInputLine(BufferedReader br) {
		String input = null;

		try {

			while (input == null) {
				System.out.println("enter input followed by enter");
				input = br.readLine();
			}

		} catch (IOException ioe) {
			System.out.println("getInputLine(): could not read from keyboard");
			return null;
		}

		return input;
	}//--getInputLine()

	// Defines

	// if the value token looks like this, there was no data on this mutation
	final String DIVISION_BY_ZERO = "#DIV/0!";

	final char NUCLEOTIDE_N = 'n';

	final double ONE_THIRD = 1.0 / 3.0;

	final double VALUE_NO_MATCH = -1.0; // this value means there was no match

	final int UNIQUE_NUCLEOTIDE_COUNT = 4; // total number of unique nucleotides
										   // (a,c,g,t)

	final char[] UNIQUE_NUCLEOTIDE_ARRAY = { 'a', 'c', 'g', 't' };

	KeyPair[] key_pairs = new KeyPair[192]; // should be 192 trinucl. with
											// probability pairs

	int key_pairs_index = 0; // the starting index of the key_pairs

	BufferedReader br = null; // a buffered reader

	public MutationSpectrum(String mutationSpectrumFileName) {

		try {
			// access the mutation spectrum file with a buffered reader
			br = new BufferedReader(new FileReader(mutationSpectrumFileName));

		} catch (FileNotFoundException fnfe) {
			throw new Error("MutationSpectrum():  " + fnfe.getMessage());
		}

		String data = "";
		String temp = "";

		try {

			// read all the data from file into a string
			temp = br.readLine();
			while (temp != null) {
				data += temp + " "; // add space between each line
				temp = br.readLine();
			}
		} catch (IOException ioe) {
			throw new Error("MutationSpectrum():  " + ioe.getMessage());
		}

		// tokenize all the data from file
		StringTokenizer st = new StringTokenizer(data);

		String token; // temporary holders
		double token_value;

		String key_token;
		double value_token;

		// collect all tokens
		while (st.hasMoreTokens()) // should have 384 (192 * 2) tokens
		{
			// get the Key token)
			token = st.nextToken();
			key_token = token;

			//System.out.println("key token = " + key_token);

			if (!st.hasMoreTokens()) {
				throw new Error(
						"MutationSpectrum constructor:  Key-Value-pair missing Value");
			}

			// read the value token
			token = st.nextToken();

			//System.out.println("next token = " + token);

			// test if there was division by zero, and thus not a legal
			// probability
			if (token.equalsIgnoreCase(DIVISION_BY_ZERO)) {
				value_token = ONE_THIRD; // the three mutations are all unlikely
										 // but will have the same unlikelyhood
			} else {
				try {
					// get the double value from token
					value_token = Double.parseDouble(token);

				} catch (NumberFormatException nfe) {
					throw new Error("Number Format Exception: "
							+ nfe.getMessage());
				}
			}//else

			// make a new key value pair and add it to array of key value pairs
			key_pairs[key_pairs_index] = new KeyPair(key_token, value_token);
			++key_pairs_index; // go to next position in array

		}//--while has more tokens

	}// constructor

	public double getTNProbability(String trinucleotide, char mutationResult) {
		// since it's a trinucleotide, length must be three (tri)
		if (trinucleotide.length() != 3) {
			throw new Error(
					"getTNProbability(): trinucleotide length not equal to THREE");
		}

		// test if the middle nucl. mutates into itself (which in this case is
		// illegal)
		if (trinucleotide.charAt(1) == mutationResult) {
			System.out.println("Error: Middle nucleotide has silent mutation");
			return VALUE_NO_MATCH;
		}

		double resultProb = VALUE_NO_MATCH;

		// test if the first nucl. is 'n', which means this trinucleotide is at
		// the very start of a sequence
		// or if the last nucl. is 'n', (at the of the sequence)
		if (trinucleotide.charAt(0) == NUCLEOTIDE_N
				|| trinucleotide.charAt(2) == NUCLEOTIDE_N) {
			if (trinucleotide.charAt(1) == NUCLEOTIDE_N)
				throw new Error("Middle trinucleotide is \"n\", illegal");

			resultProb = findEndTriNuclProb(trinucleotide, mutationResult);
			// System.out.println("average result found for end trinucleotide =
			// " + resultProb);
		} else {
			resultProb = findTriNuclProb(trinucleotide, mutationResult);

			// test for ZERO probability
			//	if(resultProb == 0.0)
			//		throw new Error("probability of lookup was ZERO");

			// System.out.println("result found for NON-end trinucleotide = " +
			// resultProb);
		}

		if (resultProb == VALUE_NO_MATCH) {
			// no match was found
			System.out.println("NO match could be found (after lookup)");

			// return value indicating no match could be found
			return VALUE_NO_MATCH;
		} else
			return resultProb;
	}//--getTNProbability()

	////////////////////////////////////////////////////////////////

	private double findTriNuclProb(String trinucleotide, char mutationResult) {
		// create the lookup key
		String lookup_key = getLookup_Key(trinucleotide, mutationResult);

		// now lookup this "new" trinucleotide
		double probability = getMatch(lookup_key);

		// return this probability
		return probability;

	}//--findTriNuclProb

	/*
	 * for the beginning and end of a sequence, there is no trinucleotide, (nxx,
	 * xxn), so we have to replace "n" with "a,c,g,t", and calculate the average
	 * probability for all four replacements
	 */
	private double findEndTriNuclProb(String generalTrinucleotide,
			char mutationResult) {
		//		if(generalTrinucleotide.charAt(0) == NUCLEOTIDE_N)
		//			System.out.println("first trinucleotide in a sequence found");
		//		if(generalTrinucleotide.charAt(2) == NUCLEOTIDE_N)
		//			System.out.println("last trinucleotide in a sequence found");
		// for all variations of n (acgt) lookup mutation and find average
		// probability

		char unique_nucl;
		double totalProbability = 0.0; // the sum of the trinucleotide
									   // probabilities
		String uniqueTrinucleotide;
		String lookup_key;

		for (int i = 0; i < UNIQUE_NUCLEOTIDE_COUNT; i++) {
			// get the unique nucleotides from array
			unique_nucl = UNIQUE_NUCLEOTIDE_ARRAY[i];

			// make a new trinucleotide by replacing the 'n' with (a,c,g,t)
			uniqueTrinucleotide = generalTrinucleotide
					.replace('n', unique_nucl);

			// create the lookup key
			lookup_key = getLookup_Key(uniqueTrinucleotide, mutationResult);

			// now lookup this "new" trinucleotide
			double curr_prob = getMatch(lookup_key);

			// sum the probs
			totalProbability += curr_prob;
		}

		// now calculate the average probability of the lookup
		double averageProbability = totalProbability / UNIQUE_NUCLEOTIDE_COUNT;

		return averageProbability;

	}//--findFirstTriNuclProb()

	////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////

	private String getLookup_Key(String trinucleotide, char mutationResult) {

		// turn the trinucl. and the mutation nucl. into a lookup key (string)
		String lookup_key = trinucleotide.charAt(0) + "("
				+ trinucleotide.charAt(1) + "->" + mutationResult + ")"
				+ trinucleotide.charAt(2);

		return lookup_key;

	}//--getLookup_Key()

	///////////////////////////////////////////////////////////////////

	/*
	 * search through all key pairs for a match with variable "lookup_key"
	 * returns probability of match if there is one else returns final value
	 * "VALUE_NO_MATCH"
	 */
	private double getMatch(String lookup_key) {
		// search through every key pair for a match
		KeyPair keyPair;
		for (int i = 0; i < key_pairs.length; i++) {
			keyPair = key_pairs[i];
			if (keyPair.key.equalsIgnoreCase(lookup_key)) {
				// there is match, so return the probability value for this
				// mutation
				//			System.out.println("value found: " + keyPair.value);
				return keyPair.value;
			}
		}//--for every key pair

		// no match found, so return value indicating no match could be found
		return VALUE_NO_MATCH;

	}//--getMatch()

	////////////////////////////////////////////////////////////////////////

	// a key pair class holding a trinucleotide mutation and it's probability
	class KeyPair {
		String key; // the trinucleotide and it's mutation as a String

		double value; // the probability of this tri nucl. mutation

		// constructor
		public KeyPair(String key, double value) {
			this.key = key;
			this.value = value;
		}

	}//--class KeyPair

	///////////////////////////////////////////////////////////////////

}//--class

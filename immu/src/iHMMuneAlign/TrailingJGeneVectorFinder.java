package iHMMuneAlign;

/**
 * rewrite By Jessica Ye
 * @author harris
 */

import java.util.*;
import org.biojava.bio.seq.*;
import gnu.bioinformatics.jaligner.JAligner;
import gnu.bioinformatics.jaligner.util.Alignment;
import java.io.*;

/**
 * contains methods for identifying the JGene and remove the trailing C or
 * Vector Region
 */
public class TrailingJGeneVectorFinder {
	/**
	 * main method for testing only
	 */
	public static void main(String[] args) {
		TrailingJGeneVectorFinder TrailingJRegionFinder = new TrailingJGeneVectorFinder();

		// minimum length of of UMS sequence after trailing J region is removed
		final int MIN_UMS_ALIGNMENT_START_OFFSET = 70;

		String mySequence = "CAGGTGCAGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGCCCCTCACCTGCGGTGTCTATGGTGGGTCCTTCACTGGTGACTTCTGGACCTGGATCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGGAAATCTATCAAAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAATAGCCACGTCCAAGAACCAATTCTCCCTGAGGCTGAATTCTTTGACCGCCGCGGACACGGCCAAATATTTCTGTGCGAGAGGCCTCTCGAATACTGCAGGTCGTCGGGGCCCACCCGCTAAGGCTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA";
		String result = TrailingJRegionFinder
				.getResult( mySequence,				
						MIN_UMS_ALIGNMENT_START_OFFSET);
		System.out.println("Best Jgene String: "+TrailingJRegionFinder.bestJGeneString);

	}//--main


	final String DEFAULT_JGENE_FILE_NAME = "/srvr/ihmmune/tomcat7/webapps/immuFile/IGHJ_Repertoire.fa";

	// gap penalties for use by JAligner
	final float OPEN_GAP_COST = (float) 10.0;
	final float EXTEND_GAP_COST = (float) 0.5;

	// score matrix to be used by the JAligner
	final String SCORE_MATRIX_NAME = "MATCH";

	// how should program behave when gap(s) are found
	final int GAPS_ALLOWED = 1; 
	final int GAPS_ERROR = 2;
	final int GAPS_IGNORE_SEQUENCE = 3; 
    /* Do we allow gap in V and J align or not? one said gap allow, the other say not allow */
	int Gap_Behaviour = GAPS_ALLOWED; 

	Sequence[] JGenes;
	String bestJGeneString;
	String bestJGeneName;

	/**
	 * constructor, a default fasta format file holding the JGenes is used
	 */
	public TrailingJGeneVectorFinder() {
		// load the JGenes from file
		Sequence[] JGenes = loadJGenes(DEFAULT_JGENE_FILE_NAME);
		this.JGenes = JGenes;

	}//--constructor

	/**
	 * finds the best matching JGene and returns the UMS string with the
	 * trailing JGene Vector/C-region part removed or returns null if no JGene
	 * could be identified or if no trailing C_region or Vector was identified
	 */
	public String getResult(String UMS_sequence,
			int MIN_UMS_ALIGNMENT_START_OFFSET) {
		
		JAligner jaligner = new JAligner(); // Smith Waterman aligner
		String UMS_string_uppercase = UMS_sequence.toUpperCase(); 
		final int UMS_LENGTH = UMS_string_uppercase.length();
		Sequence curr_JGene; // the current JGene in a Sequence format
		String curr_JGene_string; // the current JGene in a String format
		String best_alignment_JGene = null; // the whole best matching JGene
    	String best_JGene_name = null;
		Alignment curr_alignment = null;		
		Alignment best_scoring_alignment = null; // the alignment which had the best score
		float best_score = Float.MIN_VALUE; // the higher score the better, so init to lowest possible value
		float curr_score; // score found in current alignment
		int curr_alignment_gaps; // number of gaps found in current alignment

		// for every JGene
		for (int i = 0; i < JGenes.length; i++) {
	
			curr_JGene = JGenes[i];
			curr_JGene_string = curr_JGene.seqString();
			// perform alignment and get the resulting alignment object
			try {
				curr_alignment = getAlignment(jaligner, UMS_string_uppercase,
						curr_JGene_string);
			} catch (Exception ex) {
				throw new Error("do alignment failed " + ex.getMessage());
			}

			// look at the number of gaps, and the Gap behaviour, before we
			// decide if this sequence is a candidate for the best alignment
			curr_alignment_gaps = curr_alignment.getGaps(); // number of gaps
															// from current
															// alignment
			boolean sequence_accepted = true;

			if (curr_alignment_gaps > 0) // we have gaps
			{
				// decide what action to take based the set "Gap_Behaviour"
				switch (Gap_Behaviour) {
				case GAPS_ALLOWED:
					sequence_accepted = true;
					break;
				case GAPS_IGNORE_SEQUENCE:
					sequence_accepted = false;
					break;
				case GAPS_ERROR:
					throw new Error(
							"Gaps found in sequence, when Gaps Behaviour equals GAPS_ERROR");
				}//--switch
			}//--if gaps

			curr_score = curr_alignment.getScore(); // get current score from
													// alignment result

			if (curr_score > best_score && sequence_accepted) {		
				best_JGene_name = curr_JGene.getName();
				best_alignment_JGene = curr_JGene_string;
				best_score = curr_score; // current alignment score
				best_scoring_alignment = curr_alignment;
			}//--if new best alignment score

		}//--for every Jgene

		// The best alignment has now been determined
		// number of gaps found in best alignment
		int best_alignment_gaps = best_scoring_alignment.getGaps();
		
		bestJGeneString=best_alignment_JGene;
		bestJGeneName=best_JGene_name;

		// no gaps allowed in the best alignment at this stage, because the VDJ
		// model is not
		// built to accept gaps at this stage
		if (best_alignment_gaps > 0)
			throw new Error("best JGene alignment has gaps");

		// number of matching nucleotides in best alignment
		int best_alignment_similarity = best_scoring_alignment.getSimilarity();

		// the starting row and column offset found in the best alignment (row
		// is UMS sequence and column is JGene sequence)
		// its row and column becuase its a matrix
		int best_alignment_row_offset = best_scoring_alignment.getRowOffset();
		int best_alignment_col_offset = best_scoring_alignment.getColOffset();

		// find the length of the matching sequence part of the JGene
		int matchingJGenePartLength = best_alignment_JGene.length()
				- best_alignment_col_offset;

		// find the position in the UMS where the matcing part of the J Gene
		// ends, and the trailing Vector starts
		int endJGeneMatchPosition = best_alignment_row_offset
				+ matchingJGenePartLength;
		
	
		// display the the best results
		System.out.println();
		System.out.println("Best JGene match found: " + best_JGene_name);
		float score = best_scoring_alignment.getScore();
		System.out.println("score = " + score);
		System.out.println("row offset (UMS) = " + best_alignment_row_offset);
		System.out.println("col offset (JGene) = " + best_alignment_col_offset);
		System.out.println("JGene Length = " + best_alignment_JGene.length());
		

		// if the JGene ends before the end of the UMS sequence, then
		// remove the "trailing" part of the UMS string (the part after the end
		// of the JGene), (the Vector or C-Region)		
		String UMSminusTrailingJGenePart;
		if(endJGeneMatchPosition < UMS_LENGTH)
		{
			UMSminusTrailingJGenePart = UMS_sequence.substring(0,
				endJGeneMatchPosition);
			System.out.println("Trailing JGene Vector or C-region removal");
		}
		else
		{
			UMSminusTrailingJGenePart = UMS_sequence; // no removal, because the entire JGene is not present in the UMS sequence
			System.out.println("NO Trailing JGene Vector or C-region removal");
			return null; // signifies that there has is no C-Region or Vector trailing the JGene, and that part of the JGene may be missing
		}

		System.out.println("UMS minus trailing JGene part = \n"
				+ UMSminusTrailingJGenePart);
		
		// make the UMS string lowercase in the PostAlignmentResult (becuase the
		// J and D genes will be lowercase as well)
		String UMS_string_lowercase = UMS_string_uppercase.toLowerCase();

		// return the UMS string with the trailing Vector or C-region removed
		return UMSminusTrailingJGenePart;

	}//--getResult()

	/**
	 * this method calls "doAlignment()" to perfrom the actual alignment, but if
	 * the call fails (throws an exception), the method "doAlignmentReverse" is
	 * called instead.
	 * 
	 * returns the alignment object
	 */
	private Alignment getAlignment(JAligner jaligner,
			String UMS_string_uppercase, String JGene_string) {
		// get input from keyboard for pausing
		BufferedReader streami = new BufferedReader(new InputStreamReader(
				System.in));

		// the result of whatever alignment works (straight or reverse)
		Alignment result = null;
		try {
			// perform alignment and get the resulting alignment object
			result = doAlignment(jaligner, UMS_string_uppercase, JGene_string);
		} catch (java.lang.ArrayIndexOutOfBoundsException oobex) {

			System.out
					.println("Array was out of bounds: " + oobex.getMessage());

			// the regular way of aligning (UMS as row, and JGene as column
			// failed, so we'll try the reverse)
			try {
				System.out
						.println("doAlignment call failed (out of bonds exception), so we try the reverse call");
				//	streami.readLine();
				result = doAlignmentReverse(jaligner, UMS_string_uppercase,
						JGene_string);

				// since the JGene and UMS array are swapped, we need to change
				// them back
				// including their names and offsets

				// get both sequences
				char[] JGene = result.getSequence1();
				char[] UMS = result.getSequence2();

				int col_offset = result.getColOffset();
				int row_offset = result.getRowOffset();

				String JGene_name = result.getName1();
				String UMS_name = result.getName2();

				// swap them
				result.setSequence1(UMS);
				result.setSequence2(JGene);

				result.setRowOffset(col_offset);
				result.setColOffset(row_offset);

				result.setName1(UMS_name);
				result.setName2(JGene_name);

			} catch (Exception ex) // the reverse alignment call failed as well,
								   // so nothing more to do
			{
				ex.printStackTrace();
				throw new Error("Jaligner exception (Second alignment call): "
						+ ex.getMessage());
			}

		}//--catch(out of bounds exception

		catch (Exception ex) 
		{
			ex.printStackTrace();
			throw new Error("Jaligner exception: " + ex.getMessage());
		}

		return result;
	}

	/**
	 * perform an alignment with "jaligner" on UMS sequence "UMS_sequence" and
	 * JGene Sequence "JGene_sequence"
	 * 
	 * returns the alignment result object
	 */
	private Alignment doAlignment(JAligner jaligner,
			String UMS_string_uppercase, String JGene_string) throws Exception {
		
		String JGene_string_uppercase = JGene_string.toUpperCase();
		Alignment alignment_result = JAligner.sw(UMS_string_uppercase,
				JGene_string_uppercase, SCORE_MATRIX_NAME, OPEN_GAP_COST,
				EXTEND_GAP_COST);
		
		return alignment_result;
	}//--doAlignment()

	/**
	 * perform an alignment with "jaligner" on JGene Sequence "JGene_sequence"
	 * and UMS sequence "UMS_sequence" in that particular order (as opposed to
	 * the UMS then JGene first)
	 * 
	 * returns the alignment result object
	 */
	private Alignment doAlignmentReverse(JAligner jaligner,
			String UMS_string_uppercase, String JGene_string) throws Exception {
		
		String JGene_string_uppercase = JGene_string.toUpperCase();
		Alignment reverse_alignment_result = JAligner.sw(
				JGene_string_uppercase, UMS_string_uppercase,
				SCORE_MATRIX_NAME, OPEN_GAP_COST, EXTEND_GAP_COST);
		
		return reverse_alignment_result;
	}//--doAlignmentReverse()

	/**
	 * read all the JGene's from file "fastaFileName" into an array list
	 * 
	 * returns Array of JGene sequences of type "Sequence"
	 */
	private Sequence[] loadJGenes(String fastaFileName) {

		ArrayList JGenes = new ArrayList(300);
		try {
			FastaReader fr = new FastaReader();
			JGenes = fr.readFile(new java.io.File(fastaFileName));
		} catch (Exception e) {
			e.printStackTrace();
			throw new Error(e.getMessage());
		}

		Sequence[] JGene_array = new Sequence[JGenes.size()];
		Sequence temp_sequence;

		// for every Sequence in arraylist JGenes
		for (int i = 0; i < JGenes.size(); i++) {
			temp_sequence = (Sequence) JGenes.get(i);
			JGene_array[i] = temp_sequence;
		}//--for(i)

		return JGene_array;

	}//--loadJGenes()

}//--BestJGeneFinder

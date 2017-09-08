package iHMMuneAlign;

import gnu.bioinformatics.jaligner.*;
import gnu.bioinformatics.jaligner.util.Alignment;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import org.biojava.bio.seq.Sequence;
import java.util.*;
import java.io.File;
/**
 * contains methods for finding the best VGene match to a CDR sequence
 */
public class BestVGeneFinder {
	/**
	 * main method for testing only
	 */
	public static void main(String[] args) {
		BestVGeneFinder VGene_finder = new BestVGeneFinder(new File("/srvr/ihmmune/tomcat7/webapps/immuFile/IGHV_Repertoire.fa"));
		final int MIN_UMS_ALIGNMENT_END_OFFSET = 10;

		PostAlignmentResult result = VGene_finder.getResult(
						 // "CAGGTGCAGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGCCCCTCACCTGCGGTGTCTATGGTGGGTCCTTCACTGGTGACTTCTGGACCTGGATCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGGGAAATCTATCAAAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAATAGCCACGTCCAAGAACCAATTCTCCCTGAGGCTGAATTCTTTGACCGCCGCGGACACGGCCAAATATTTCTGTGCGAGAGGCCTCTCGAATACTGCAGGTCGTCGGGGCCCACCCGCTAAGGCTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA",
						   "CAGGTGCAGCTAAAACAGTGGGGCGCAGGACTGNTGAAGCCTTCGGAGGCCCTGTCCCACACCTGCGGTGTCTATGGTGGGTCCTTCTCTGGTTACTTCTGGACCTGGATCCGCCAGGTCCCAGGGAGGGGGCTGGAGTGGATTGGGGAAATCAATCAAAGTGGAAGCACCAAGTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAATAGACACGTCTAAGAACCACTTCTCCCTGCGGCTGAGTTCAATGACCGCCGCGGACACAGCTGAATATTTCTGTGCGAGAGGCCTTCCGGGTACTGCAGGTCGTCGGGGCCCACCCGCTAAGGCTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA",							
								MIN_UMS_ALIGNMENT_END_OFFSET);
         result.displayAlignmentInfo();
	}//--main

	
	final float OPEN_GAP_COST = (float) 10.0;
	final float EXTEND_GAP_COST = (float) 0.5;
	final String SCORE_MATRIX_NAME = "MATCH";

	final int GAPS_ALLOWED = 1; 
	final int GAPS_ERROR = 2; 
	final int GAPS_IGNORE_SEQUENCE = 3; 	
	int Gap_Behaviour = GAPS_ALLOWED;  //allow gaps_in_alignment

	// array of VGenes for alignment with UMS sequence
	Sequence[] VGenes;

	/**
	 * constructor, "fileName" is the name of a file holding all the fasta
	 * format VGene sequences
	 */
	public BestVGeneFinder(File fastaFile) {
		Sequence[] VGenes = loadVGenes(fastaFile);
		this.VGenes = VGenes;

	}//--constructor

	/**
	 * finds the best matching VGene and return a "PostAlignmentResult" with the
	 * original UMS sequence and best matching VGene sequence. There must be at
	 * least "MIN_UMS_ALIGNMENT_END_OFFSET" number of nucleotides trailing the
	 * alignment part of the UMS string, because this is where the D and J gene
	 * will sit
	 */
	public PostAlignmentResult getResult(String UMS_sequence,
			int MIN_UMS_ALIGNMENT_END_OFFSET) {
		
		JAligner jaligner = new JAligner(); 

		String UMS_string_uppercase = UMS_sequence.toUpperCase(); 
		final int UMS_LENGTH = UMS_string_uppercase.length(); 

		Sequence curr_VGene; 
		String curr_VGene_string; 
		String best_alignment_VGene = null; 
		String best_VGene_name = null;
		Alignment curr_alignment = null;
		Alignment best_scoring_alignment = null;

		float best_score = Float.MIN_VALUE; 
		float curr_score; 
		int curr_alignment_gaps; // number of gaps found in current alignment

		// for every VGene
		for (int i = 0; i < VGenes.length; i++) {
			curr_VGene = VGenes[i];
			curr_VGene_string = curr_VGene.seqString();

			// perform alignment and get the resulting alignment object
			try {
				curr_alignment = getAlignment(jaligner, UMS_string_uppercase,
						curr_VGene_string);
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
					throw new Error("Gaps found in sequence, when Gaps Behaviour equals GAPS_ERROR");
				}//--switch
			}//--if gaps

			// check if the aligned part (from alignment til end) of the UMS is
			// at least X nucl. longer than
			// the aligned part (from alignment til end) of the VGene
			int from_alignment_length_UMS = UMS_LENGTH
					- curr_alignment.getRowOffset();
			int curr_VGene_length = curr_VGene_string.length()
					- curr_alignment.getColOffset();
			int from_alignment_length_VGene = curr_VGene_length
					- curr_alignment.getRowOffset();

			if (from_alignment_length_VGene > (from_alignment_length_UMS - MIN_UMS_ALIGNMENT_END_OFFSET)) {
				sequence_accepted = false;
			}//--if

			curr_score = curr_alignment.getScore(); // get current score from
													// alignment result

			if (curr_score > best_score && sequence_accepted) { 			
				try {
					BufferedReader streami = new BufferedReader(
							new InputStreamReader(System.in));
				
				} catch (Exception ex) {
				}

				best_VGene_name = curr_VGene.getName();
				best_alignment_VGene = curr_VGene_string;
				best_score = curr_score; 
				best_scoring_alignment = curr_alignment;

			}//--if new best alignment score

		}//--for every Vgene

	
		int best_alignment_gaps = best_scoring_alignment.getGaps();

		// no gaps allowed in the best alignment at this stage, because the VDJ
		// model is not built to accept gaps at this stage
		if (best_alignment_gaps > 0)
			return null; // null means that the best alignment has gaps

		// number of matching nucleotides in best alignment
		int best_alignment_similarity = best_scoring_alignment.getSimilarity();

		// the row and column offset found in the best alignment (row is UMS
		// sequence and column is VGene sequence)
		int best_alignment_row_offset = best_scoring_alignment.getRowOffset();
		int best_alignment_col_offset = best_scoring_alignment.getColOffset();

		// get the names of the UMS sequence and the VGene sequence
		
		String UMS_name = "not noted";
		String UMS_string_lowercase = UMS_string_uppercase.toLowerCase();

		// make a postalignment object which holds the most important
		// information about the alignment
		PostAlignmentResult result = new PostAlignmentResult(
				UMS_string_lowercase, best_alignment_VGene,
				best_alignment_row_offset, best_alignment_col_offset,
				best_alignment_gaps, best_alignment_similarity, UMS_name,
				best_VGene_name);

		// display the best results
		float score = best_scoring_alignment.getScore();
		System.out.println("score = " + score);
		System.out.println("row offset = " + best_alignment_row_offset);
		System.out.println("col offset = " + best_alignment_col_offset);
		
		// return this object
		return result;

	}//--getResult()

	/**
	 * this method calls "doAlignment()" to perfrom the actual alignment, but
	 * the call fails (throws an exception), the method "doAlignmentReverse" is
	 * called instead.
	 * 
	 * returns the alignment object
	 */
	private Alignment getAlignment(JAligner jaligner,
			String UMS_string_uppercase, String VGene_string) {
		// get input from keyboard for pausing
		BufferedReader streami = new BufferedReader(new InputStreamReader(
				System.in));

		// the result of whatever alignment works (straight or reverse)
		Alignment result = null;
		try {
			// perform alignment and get the resulting alignment object
			result = doAlignment(jaligner, UMS_string_uppercase, VGene_string);
		} catch (java.lang.ArrayIndexOutOfBoundsException oobex) {
			System.out
					.println("Array was out of bounds: " + oobex.getMessage());

			// the regular way of aligning (UMS as row, and VGene as column
			// failed, so we'll try the reverse)
			try {
				System.out
						.println("doAlignment call failed (out of bonds exception), so we try the reverse call");
				result = doAlignmentReverse(jaligner, UMS_string_uppercase,
						VGene_string);

				char[] VGene = result.getSequence1();
				char[] UMS = result.getSequence2();

				int col_offset = result.getColOffset();
				int row_offset = result.getRowOffset();

				String VGene_name = result.getName1();
				String UMS_name = result.getName2();

				// swap them
				result.setSequence1(UMS);
				result.setSequence2(VGene);

				result.setRowOffset(col_offset);
				result.setColOffset(row_offset);

				result.setName1(UMS_name);
				result.setName2(VGene_name);

			} catch (Exception ex) 
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

		// return the alignment result "result"
		return result;
	}

	/**
	 * perform an alignment with "jaligner" on UMS sequence "UMS_sequence" and
	 * VGene Sequence "VGene_sequence"
	 * 
	 * returns the alignment result object
	 */
	private Alignment doAlignment(JAligner jaligner,
			String UMS_string_uppercase, String VGene_string) throws Exception {	
			
		String VGene_string_uppercase = VGene_string.toUpperCase();
		Alignment alignment_result = JAligner.sw(UMS_string_uppercase,
				VGene_string_uppercase, SCORE_MATRIX_NAME, OPEN_GAP_COST,
				EXTEND_GAP_COST);
		
		return alignment_result;
	}//--doAlignment()

	/**
	 * perform an alignment with "jaligner" on VGene Sequence "VGene_sequence"
	 * and UMS sequence "UMS_sequence" in that particular order (as opposed to
	 * the UMS then VGene first)
	 * 
	 * returns the alignment result object
	 */
	private Alignment doAlignmentReverse(JAligner jaligner,
			String UMS_string_uppercase, String VGene_string) throws Exception {
		
		String VGene_string_uppercase = VGene_string.toUpperCase();
		Alignment reverse_alignment_result = JAligner.sw(
				VGene_string_uppercase, UMS_string_uppercase,
				SCORE_MATRIX_NAME, OPEN_GAP_COST, EXTEND_GAP_COST);
		
		return reverse_alignment_result;
	}//--doAlignment()

	/**
	 * read all the VGene's from file "fastaFileName" into an array list
	 * 
	 * returns Array of VGene sequences of type "Sequence"
	 */
	private Sequence[] loadVGenes(java.io.File fastaFile) {

		ArrayList VGenes = new ArrayList(300);
		try {
			FastaReader fr = new FastaReader();
			VGenes = fr.readFile(fastaFile);
		} catch (Exception e) {
			
			e.printStackTrace();
			throw new Error(e.getMessage());
		}

		Sequence[] VGene_array = new Sequence[VGenes.size()];
		Sequence temp_sequence;

		// for every Sequence in arraylist VGenes
		for (int i = 0; i < VGenes.size(); i++) {
			temp_sequence = (Sequence) VGenes.get(i);
			VGene_array[i] = temp_sequence;
		}//--for(i)

		return VGene_array;

	}//--loadVGenes()

}//--BestVGeneFinder

/* this class used to retrieve results, can add to Jgene */
/**
 * this class holds the amputated (start offset removed) UMS string and VGene
 * string
 */

class PostAlignmentResult {
	
	private String UMS_string;
	private String VGene_string;
	private String UMS_name;
	private String VGene_name;
	private int row_offset;
	private int column_offset;
	private int gaps_in_alignment;
	private int similarity;

	/**
	 * constructor taking original unmodified UMS and VGene string, the row
	 * offset relates to the UMS_string, and the column offset relates to the
	 * VGene_string
	 */
	public PostAlignmentResult(String UMS_string, String VGene_string,
			int row_offset, int column_offset, int gaps_in_alignment,
			int similarity, String UMS_name, String VGene_name) {
		this.UMS_string = UMS_string;
		this.VGene_string = VGene_string;
		this.row_offset = row_offset;
		this.column_offset = column_offset;
		this.gaps_in_alignment = gaps_in_alignment;
		this.similarity = similarity;
		this.UMS_name = UMS_name;
		this.VGene_name = VGene_name;

	}//--PostAlignmentResult

	/**
	 * get the UMS string from the start of the alignment til the end The number
	 * of characters in the row offset is removed from the start of this string
	 */
	public String getUMSfromAlignmentString() {
			String UMSfromAlignment = UMS_string.substring(row_offset, UMS_string
				.length());

		return UMSfromAlignment;

	}//--getUMSfromAlignmentString()

	/**
	 * get the VGene string from the start of the alignment til then end The
	 * number of characters in the column offset is removed from the start of
	 * this string
	 */
	public String getVGeneAlignmentString() {
		String VGeneFromAlignment = VGene_string.substring(column_offset,
				VGene_string.length());

		return VGeneFromAlignment;

	}//--getVGeneAlignmentString()

	/**
	 * get method for name UMS sequence
	 */
	public String getUMSname() {
		return UMS_name;
	}//--getUMSname()

	/**
	 * get method for name VGene sequence
	 */
	public String getVGeneName() {
		return VGene_name;
	}//--getUMSname()

	/**
	 * get method for the column offset
	 */
	public int getColOffset() {
		return this.column_offset;
	}//--getColOffset()

	/**
	 * get method for the column offset
	 */
	public int getRowOffset() {
		return this.row_offset;
	}//--getRowOffset()

	/**
	 * display information about this alignment
	 */
	public void displayAlignmentInfo() {
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println("*** VGene ALIGNMENT INFO ***");
		System.out.println();

		// display info on UMS string

		//System.out.println("UMS name: " + this.UMS_name);
		System.out.println("complete UMS string length = "
				+ UMS_string.length());
		System.out.println("complete UMS string:");
		System.out.println(UMS_string);

		System.out.println();

		String from_alignment_UMS = this.getUMSfromAlignmentString();
		System.out.println("from alignment UMS length: "
				+ from_alignment_UMS.length());
		System.out.println("UMS offset: " + row_offset);
		System.out.println("from alignment UMS:");
		System.out.println(from_alignment_UMS);
		System.out.println();

		// display info on VGene string

		System.out.println("VGene: " + this.VGene_name);
		String from_alignment_VGene = this.getVGeneAlignmentString();
		System.out.println("from alignment VGene length: "
				+ from_alignment_VGene.length());
		System.out.println("VGene offset: " + column_offset);
		System.out.println("from alignment VGene:");
		System.out.println(from_alignment_VGene);

		// display number of alignments found (nucl. matches)
		System.out.println();
		System.out.println("similarity: " + similarity);

		// display number of gaps found in alignment
		System.out.println("Gaps in alignment: " + gaps_in_alignment);

		System.out.println();
		System.out.println();

	}//--displayAlignmentInfo()

}//--PostAlignmentResult

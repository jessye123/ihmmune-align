package iHMMuneAlign;


import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPFactory;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.StatePath;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;

import java.io.*;
import immu.IHMMuneResult;


public class AlignmentThread extends Thread implements GlobalDefines {
   
	IHMMuneResult myIHMMue; //this hold the iHHMune Align result
	String AlignMatchResult;
	
	Sequence sequence;

	ProbabilityHolder probHolder;

	MutabilityScore mutabilityScore;

	File V_GENES_FILE, D_GENES_FILE, J_GENES_FILE;

	byte alignmentType;

	int dGeneAcceptanceType;
	
	PrintWriter writer;
	
	

	public AlignmentThread(PrintWriter writer, Sequence sequence, File V_GENES_FILE, File D_GENES_FILE, File J_GENES_FILE,
			ProbabilityHolder probHolder, MutabilityScore mutabilityScore,
			 byte alignmentType, int dGeneAcceptanceType) {
		this.writer =writer;
		this.sequence = sequence;
		this.probHolder = probHolder;
		this.mutabilityScore = mutabilityScore;
		this.V_GENES_FILE = V_GENES_FILE;
		this.D_GENES_FILE = D_GENES_FILE;
		this.J_GENES_FILE = J_GENES_FILE;
		this.alignmentType = alignmentType;
		this.dGeneAcceptanceType = dGeneAcceptanceType;
		
	}

	
    //this is important method //Test at MultiCellTool.java and SingleCellTool.java
	public Object runAlignment() {

		String sequenceString = removeFastaStyleWhiteSpace(sequence.seqString());
		String sequenceName = (String) sequence.getName();
		//write to file
		if(writer!=null){
				writer.println();
				writer.print(sequenceString);  
			}
		

// ***************************First identify the VGene in this Sequence*********************************
		
		BestVGeneFinder vfinder = new BestVGeneFinder(V_GENES_FILE);
		
		final int MIN_UMS_ALIGNMENT_END_OFFSET = 12; //VGene length must be  "MIN_UMS_ALIGNMENT_END_OFFSET" nucleotides shorter than the seq
		
		PostAlignmentResult alignment_results = vfinder.getResult( sequenceString, MIN_UMS_ALIGNMENT_END_OFFSET);

		if (alignment_results == null) {
			throw new Error("VGene finder result equal NULL");			
		}

		String from_alignment_start_UMS_string = alignment_results.getUMSfromAlignmentString();
		String from_alignment_start_VGene_string = alignment_results.getVGeneAlignmentString();
		String UMS_name = alignment_results.getUMSname();
		String VGene_name = alignment_results.getVGeneName();

		Sequence from_alignment_start_UMS_seq = null;
		Sequence from_alignment_start_VGene_seq = null;
		try {
			// create UMS and VGene "Sequence" based on sequence string and name
			from_alignment_start_UMS_seq = DNATools.createDNASequence(from_alignment_start_UMS_string, UMS_name);
			from_alignment_start_VGene_seq = DNATools.createDNASequence(from_alignment_start_VGene_string, VGene_name);
		} catch (IllegalSymbolException illse) {
			throw new Error(illse.getMessage());
		}
		
		int VGene_start_offset = alignment_results.getColOffset();		
		int completeVGeneLength = VGene_start_offset + from_alignment_start_VGene_string.length();

		/*
		 * ***************** A-value calculation *****************
		 */
		// calculate the A value of this VGene using the common region and
		// number of mutations found
		// if the VGene does not contain the entire Common Region, an Error will be thrown
		/* This A_Score will be OK to identify number of mutation in V gene, and also find CDR regions, CDR1+CDR2, How to tell apart???? */
		// however identify CDR by 128-254 maybe a mistake
		double A_probability = A_Score.A_probability(
				from_alignment_start_VGene_string,
				from_alignment_start_UMS_string, VGene_start_offset, null);

		System.out.println("A_probability = " + A_probability);
		
		
// ***************************ABOVE Vgene*******Below Jgene Jaligner find best J gene and chop C region*************		
        
				// remove the possibly trailing J gene Vector part (or C region)
				TrailingJGeneVectorFinder vectorFinder = new TrailingJGeneVectorFinder();
				String UMS_no_C = vectorFinder.getResult(from_alignment_start_UMS_string, 70);

				String from_alignment_start_UMS_no_C_string = null;
				Sequence from_alignment_start_UMS_no_C_seq = null;

				if (UMS_no_C != null) {				
					System.out.println("C_Region Removal");			
					from_alignment_start_UMS_no_C_string = UMS_no_C;
					try {
						from_alignment_start_UMS_no_C_seq = DNATools.createDNASequence(UMS_no_C, "UMS_name");
					} catch (IllegalSymbolException illse) {
						throw new Error(illse.getMessage());
					}
					
				} else {				
					System.out.println("NO C_Region Removal");
					from_alignment_start_UMS_no_C_string = from_alignment_start_UMS_string;
					from_alignment_start_UMS_no_C_seq = from_alignment_start_UMS_seq;
				}
				
				 Sequence JgeneSeq;
				try {
					JgeneSeq = DNATools.createDNASequence(vectorFinder.bestJGeneString, vectorFinder.bestJGeneName);
					
				} catch (IllegalSymbolException illse) {
					throw new Error(illse.getMessage());
				}
		//Above identify correct Jgene It will be waste to use Markov model to do this again! program eliminated repeated loops
		// ***************************ABOVE Vgene, Jgene and chop C region*******Below Markov model align*************
				
				// a reference to the class which creates the Markov Model
				VpnpDpnpJCnoC vdj;

				// a MarkovModel reference (the result of creating a Markov Model)
				MarkovModel markov_model;	
			
				vdj = new VpnpDpnpJCnoC(mutabilityScore);
				
				//To improve efficient, Changed J_GENES_FILE to Gene seq	
				markov_model = vdj.createModel(from_alignment_start_VGene_seq,
						VGene_start_offset, completeVGeneLength, D_GENES_FILE, JgeneSeq,
						probHolder, false, A_probability);

				
				// Perform Dynamic Programming algorithms on HMM
				System.out.println("about to create DP");
				DP dp = null;
				try {
					dp = DPFactory.DEFAULT.createDP(markov_model);
				} catch (BioException bioe) {
					throw new Error(bioe.getMessage());
				}
				System.out.println("Dp created");

				// perform Viterbi algorithm on Sequence/Model
				SymbolList[] res_array = { from_alignment_start_UMS_no_C_seq };
			
				System.out.println("about to create state path viterbi");
				StatePath v = null;
				try {
					v = dp.viterbi(res_array, ScoreType.PROBABILITY);
				} catch (IllegalAlphabetException illalphe) {
					throw new Error(illalphe.getMessage());
				} catch (IllegalSymbolException illsymbole) {
					throw new Error(illsymbole.getMessage());
				} catch (IllegalTransitionException illsymbole) {
					throw new Error(illsymbole.getMessage());
				}

				double statePathProbability = v.getScore();
				System.out.println("state path viterbi created");
				SymbolList viterbi_state_seq = v.symbolListForLabel(StatePath.STATES);

		//*******************display State Information**************************************

				// display information on emitting states for each nucleotide, comment out if no need
				displayStateInfo(viterbi_state_seq,
						from_alignment_start_UMS_no_C_string, sequenceName);

		//*******************display  Emission Information,write output Text file*********************
		// display information on the details of each region (Gene,P and N)
				//Old version displayEmissionInfo method
			  	
				System.out.println("\n" + "*** Details ***" + "\n");

				// display the state path probability
				System.out.println("State Path Probability = " + statePathProbability);
				

				// first, print the VGene Properties
			    String UMS = from_alignment_start_UMS_no_C_string;
				GeneInfo VGeneInfo = displayGene( UMS, viterbi_state_seq,
						V_STATE_TOKEN);
				RegionInfo p1Region = displayNRegion(UMS,
						viterbi_state_seq, V_END_P_TOKEN);
				RegionInfo n1Region = displayNRegion(UMS,
						viterbi_state_seq, VD_N_TOKEN);
				RegionInfo p2Region = displayNRegion(UMS,
						viterbi_state_seq, D_START_P_TOKEN);
				
				String n1="";
				if(p1Region!=null) n1+=p1Region.region_string;
				if(n1Region!=null) n1+=n1Region.region_string;
				if(p2Region!=null) n1+=p2Region.region_string;
				
				GeneInfo DGeneInfo = displayGene(UMS, viterbi_state_seq,
							D_STATE_TOKEN);
				
				RegionInfo p3Region = displayNRegion( UMS,
						viterbi_state_seq, D_END_P_TOKEN);
				RegionInfo n2Region = displayNRegion(UMS,
						viterbi_state_seq, DJ_N_TOKEN);
				RegionInfo p4Region = displayNRegion(UMS,
						viterbi_state_seq, J_START_P_TOKEN);

				String n2="";
				if(p3Region!=null) n2+=p3Region.region_string;
				if(n2Region!=null) n2+=n2Region.region_string;
				if(p4Region!=null) n2+=p4Region.region_string;
			
				GeneInfo JGeneInfo = displayGene(UMS, viterbi_state_seq,
						J_STATE_TOKEN);
				
				// write result to an OUTPUT file
				if(writer!=null){
				writer.print(","+VGeneInfo.gene_name); writer.print(","+DGeneInfo.gene_name);  writer.print(","+JGeneInfo.gene_name); 
				writer.print(","+VGeneInfo.aligned_ums_string); writer.print(","+n1); writer.print(","+DGeneInfo.aligned_ums_string); writer.print(","+n2); writer.print(","+JGeneInfo.aligned_ums_string); 
				writer.print(","+VGeneInfo.getMutation()); writer.print(","+DGeneInfo.getMutation());  writer.print(","+JGeneInfo.getMutation()); 
				}			
					  
				myIHMMue = new IHMMuneResult(sequenceString, VGeneInfo.gene_name,DGeneInfo.gene_name,JGeneInfo.gene_name,
						VGeneInfo.aligned_ums_string,n1,DGeneInfo.aligned_ums_string,n2,JGeneInfo.aligned_ums_string,
						Integer.toString(VGeneInfo.getMutation()),Integer.toString(DGeneInfo.getMutation()), Integer.toString(JGeneInfo.getMutation()) );
				
//*******************PostAlignment   D gene accept, Stop Codon check and Jgene in Frame *********************
				
				System.out.println("******* Begin PostAlign *******\n ");
				
				PostAlignProcessing postAlign = new PostAlignProcessing();
			
				// D gene accept
				boolean acceptedDGene = postAlign.acceptDGene(DGeneInfo, dGeneAcceptanceType);
				if (!acceptedDGene) {
					DGeneInfo.acceptedAlignment = false;
				}
				System.out.println("D gene accept : " + acceptedDGene +"\n" );
								
				
				// test if JGene is read in-frame
				JGeneReadingFrame jGeneReadingFrame = new JGeneReadingFrame(UMS, JGeneInfo, VGeneInfo);
				boolean isJGeneInFrame = jGeneReadingFrame.isJGeneInFrame();
				int relativeMotifPosition = jGeneReadingFrame.getWGXGmotifInAlignedJGenePos();
				System.out.println();
			
	  //******************* hint to Garabage collector***********************************************************
				System.gc();
				System.out.println("HINTING TO GC");
		

		// there was no Exceptions, so return null
		return myIHMMue ;// return the result 

	}//--runAlignment
    
	public IHMMuneResult getIHMMuneResult (){
		return myIHMMue;
	}
	
	public String getAlignmatchResult (){
		return AlignMatchResult;
	}
	////////////////////////////////////////////////////////////////////////

	/**
	 * display information about the alignment and what state emitted each
	 * symbol. This info is displayed in the JTextArea "outputArea"
	 */
	private void displayStateInfo(
			SymbolList viterbi_state_seq, String UMS, String selected_ums_name) {

		char UMSEmit;
		char geneToken;
		char statePreEmit;

		// strings displaying the result
		String stateResult = "";
		String umsResult = "";
		String preEmitResult = "";
		String matchResult = "";
		
		int ums_index = 0;
		Annotation annotation;
		StateInfo si;

		// ums string is shorter than the state path because of the dot states emitted
		for (int i = 1; i <= viterbi_state_seq.length(); i++) {
			annotation = viterbi_state_seq.symbolAt(i).getAnnotation();
			si = (StateInfo) annotation.getProperty(null);
			geneToken = si.geneToken;
			statePreEmit = si.preEmissionSymbol;

			// display everything except for dot states
			if (geneToken != DOT_STATE_TOKEN) {
				
				//state result
				if (geneToken == DJ_N_TOKEN || geneToken == VD_N_TOKEN) {
					stateResult += N_TOKEN;  
				} else if (geneToken == V_END_P_TOKEN
						|| geneToken == D_START_P_TOKEN
						|| geneToken == D_END_P_TOKEN
						|| geneToken == J_START_P_TOKEN) {
					stateResult += P_TOKEN;  
				} else {
					stateResult += geneToken;  				
				}
				
				//gene Emit result
				preEmitResult += statePreEmit; 
				
				// ums result
				UMSEmit = UMS.charAt(ums_index++);
				umsResult += UMSEmit;

				// Match result
				if (statePreEmit == ANY_PRE_EMISSION_TOKEN) {
					matchResult += ANY_PRE_EMISSION_TOKEN;
				} else if (statePreEmit == UMSEmit) {
					matchResult += MATCH_TOKEN;
				} else {
					matchResult += MISMATCH_TOKEN;
				}
			}
		  else {
				// it's a DOT State
				// do nothing			
			}


		}//--for every state
		
		System.out.println("\n");
		System.out.println("******** Result of aligning " + selected_ums_name + "  *******" + "\n");
		System.out.println("emitting state : " + stateResult + "\n");
		System.out.println("ums region     : " + umsResult + "\n");
		System.out.println("nucl. in gene  : " + preEmitResult + "\n");
		System.out.println("match result   : " + matchResult + "\n");
		
		AlignMatchResult = "******** Result of aligning " + selected_ums_name + "  *******" + "\n"
		                  + "emitting state : " + stateResult + "\n"
		                  +"ums region     : " + umsResult + "\n"
		                  +"nucl. in gene  : " + preEmitResult + "\n"
		                  +"match result   : " + matchResult + "\n";
		

	}

	//--displayStateInfo

	////////////////////////////////////////////////////////////////////////

	/**
	 * displays the details of an alignment (V D J genes)
	 */
	private GeneInfo displayGene(String ums,
			SymbolList viterbi_state_seq, char state_token) {
	
		Annotation anno;
		StateInfo si;
		char geneToken;

		String GeneName = "";
		String Gene = "";
	    String UMS = "";		
		String Match = "";
		int first_nucleotide_number = -1;
		int last_nucleotide_number = -1;

		int ums_index = 0;
		int state_index = 1;  // this index should match the index of the state path

	  // find the start of the Gene V/D/J State
		boolean foundStateRegion = false;
		do {
			anno = viterbi_state_seq.symbolAt(state_index).getAnnotation();
			si = (StateInfo) anno.getProperty(null);
			geneToken = si.geneToken;

			if (geneToken == state_token) {
				foundStateRegion = true;  // yes, so stop looking
				first_nucleotide_number = si.getNucleotidePosition();  // offset into Gene				
				
			} else if (geneToken == DOT_STATE_TOKEN) {
				// dot states should be ignored (they are extra), so don't increment the ums count
				state_index++;
			} else {
				
				state_index++;
				ums_index++;
			}
		} while (state_index <= viterbi_state_seq.length() && !foundStateRegion);

		if (!foundStateRegion) {
			
			return null;
		}
		
	 // in the Gene V D J State
		while (geneToken == state_token) {
		
			UMS += ums.charAt(ums_index);
			Gene += si.preEmissionSymbol;
			
			if (si.preEmissionSymbol == ums.charAt(ums_index)) {
				Match += MATCH_TOKEN;
			} else {
				Match += MISMATCH_TOKEN;
			}

			state_index++;
			ums_index++;

			// make sure we dont jump off the end of the array
			if (state_index > viterbi_state_seq.length()) {
				break;
			}

			// find the next state
			anno = viterbi_state_seq.symbolAt(state_index).getAnnotation();
			si = (StateInfo) anno.getProperty(null);
			geneToken = si.geneToken;
		}

		// find the last state - and get the name from this state
		anno = viterbi_state_seq.symbolAt(state_index - 1).getAnnotation();
		si = (StateInfo) anno.getProperty(null);

		// find the last nucleotide number
		last_nucleotide_number = si.getNucleotidePosition();
		
		GeneName = si.geneName + "  length of complete gene:  "
				+ si.completeGeneLength + "  start: " + first_nucleotide_number
				+ "  end: " + last_nucleotide_number + "\n";

		// print the Result
	System.out.println(GeneName);
	System.out.println("UMS Seq   : " + UMS +   "\n");
	System.out.println("Gene Seq  : " + Gene +  "\n");
	System.out.println("Match     : " + Match + "\n");
	System.out.println("\n");

		// //////////////// create a GeneInfo object and return it
		int completeGeneLength = si.completeGeneLength;

		GeneInfo tempGeneInfo = new GeneInfo(si.geneName, Gene,
				UMS, completeGeneLength, first_nucleotide_number,
				last_nucleotide_number);

		return tempGeneInfo;
	} //--displayGene()
	
	/**
	 * display a non gene region, can be N or P region, 
	 */
	private RegionInfo displayNRegion(String ums,
			SymbolList viterbi_state_seq, char state_token) {
	
		Annotation anno;
		StateInfo si;
		char geneToken;

		String GeneName = "";
		String UMS = "";  	// hold the nucleotide emission
	
		int ums_start_location = -1;
		int ums_end_location = -1;

		int ums_index = 0;
		int state_index = 1;

	  // find the region
		boolean foundStateRegion = false;
		do {
			anno = viterbi_state_seq.symbolAt(state_index).getAnnotation();
			si = (StateInfo) anno.getProperty(null);
			geneToken = si.geneToken;

			if (geneToken == state_token) {
				foundStateRegion = true;
				ums_start_location = ums_index;
			} else if (geneToken == DOT_STATE_TOKEN) {			
				state_index++;
			} else {
				state_index++;
				ums_index++;
			}
		} while (state_index <= viterbi_state_seq.length() && !foundStateRegion);

		if (!foundStateRegion) {
			return null;
		}
		// no, so return without no further adou

		// match up all the states
		while (geneToken == state_token) {
	
			UMS += ums.charAt(ums_index);
			state_index++;
			ums_index++;

			// make sure we dont jump off the end of the array
			if (state_index > viterbi_state_seq.length()) {
				break;
			}

			// find the next state
			anno = viterbi_state_seq.symbolAt(state_index).getAnnotation();
			si = (StateInfo) anno.getProperty(null);
			geneToken = si.geneToken;
		}

		// find the ums end pos and save it
		ums_end_location = ums_index - 1;

		// find the last state - and get the name from this state
		anno = viterbi_state_seq.symbolAt(state_index - 1).getAnnotation();
		si = (StateInfo) anno.getProperty(null);

		GeneName = si.geneName + "  ums start index:  " + ums_start_location
				+ "  ums end index:  " + ums_end_location + "\n";

		// print the Result
		System.out.println(GeneName);
		System.out.println("emitted sequence : " + UMS + "\n");
		System.out.println("\n");

		// ///////////// create the RegionInfo object and return it
		RegionInfo tempRegionInfo = new RegionInfo(UMS,
				ums_start_location, ums_end_location);

		return tempRegionInfo;
	}

	//--displayNRegion()

	////////////////////////////////////////////////////////////////////////


	/**
	 * Removes White Space characters from a Sring representation of a Sequence.
	 * The white space characters appears as '-" because of the fasta reader
	 */
	private String removeFastaStyleWhiteSpace(String string) {

		String noWhiteSpaceString = ""; 
		char curr; 
		for (int i = 0; i < string.length(); i++) {
			curr = string.charAt(i);
			if (curr != '-')
				noWhiteSpaceString = noWhiteSpaceString + curr;
		}//for

		return noWhiteSpaceString; 
	}//--removeWhiteSpace()

	
}//AlignmentThread

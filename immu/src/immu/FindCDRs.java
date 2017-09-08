package immu;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import org.biojava.bio.seq.Sequence;

import gnu.bioinformatics.jaligner.JAligner;
import gnu.bioinformatics.jaligner.util.Alignment;

import iHMMuneAlign.FastaReader;

public class FindCDRs {

	/**
	 * This class find CDR1 CDR2 and CDR3
	 * IMGT Scientific Chart Amino Acid number CDR1 27-38   CDR2 56- 65  CDR3 105-116 
	 */
	// use to return the information to website
	String CDR1, CDR2, CDR3;
	String CDR1AA, CDR2AA, CDR3AA;
	String CDRSdisplay;
	String AminoAcid, AminoAcidDisplay;
	
    //main for test ONLY
	public static void main(String[] args) throws Exception {
		
	    String inputFilename= "src/output.csv";	
		File iHMMuneAlignFile = new File(inputFilename);	
		   BufferedReader br = new BufferedReader(new FileReader(iHMMuneAlignFile));
		   String delimiter = ",";

			// work through each line in the input file
		   for (String line = br.readLine(); line != null; line = br
					.readLine()) {
			if(!line.startsWith("#")){	
				line = line.replace("\"", "");
				String[] lineData = line.split(delimiter);
				
				if (lineData.length > 1) { // skip blank lines
					String seqID = lineData[0];	
					//test only
					     seqID=seqID.substring(2);	
					
					String IGHVName = lineData[1];			
					String IGHDName = lineData[2];
					String IGHJName = lineData[3];
					
					String IGHVSeq = lineData[4];
					//change for testing
					      IGHVSeq= IGHVSeq.substring(2);
					
					String N1Seq = lineData[5];
					String IGHDSeq = lineData[6];
					String N2Seq = lineData[7];
					String IGHJSeq = lineData[8];
					
					String VMut = lineData[9];
					String DMut = lineData[10];
					String JMut = lineData[11];
					IHMMuneResult myResult = new IHMMuneResult(seqID,IGHVName,IGHDName,IGHJName,IGHVSeq,N1Seq, IGHDSeq, N2Seq,IGHJSeq,VMut,DMut,JMut);	
					
					FindCDRs findcdr = new FindCDRs();
					findcdr.GetCDRs(myResult);
				}// end if
				}//end if
			}//end for

	}	
	
	public void GetCDRs(IHMMuneResult myResult){
		try { 
		String IGHVName= myResult.getIGHVName();
		String IGHVSeq= myResult.getVSeq();
		String seqID= myResult.getSeqID();
		String VgeneSeq = FindVgeneSeq(IGHVName);
		
			
		JAligner jaligner = new JAligner(); 
		Alignment Align = doAlignment(jaligner, IGHVSeq.toUpperCase(), VgeneSeq.toUpperCase());
		int Seq_offset = Align.getRowOffset();
		int Vgene_offset = Align.getColOffset();
		int VStart= Vgene_offset  -  Seq_offset;  // V gene may longer than the seq# string
		System.out.println("v gene offset minus Seq offset: " + VStart);
	
		CDR1 = getCDR1seq(IGHVSeq, VStart);		
	    CDR2 = getCDR2seq(IGHVSeq, VStart);
	    CDRSdisplay =  "CDR1 Seq          :" + CDR1+"\n"; 
	    
	    //compare with unMuatated CDR1
	    String unMutatedCDR1 = getCDR1seq(VgeneSeq, 0);		
	    String unMutatedCDR2 = getCDR2seq(VgeneSeq, 0);    	
		CodonTranslate CodonTrans = new CodonTranslate();	
		CDR1AA= CodonTrans.findCodons(CDR1, 0);
		CDRSdisplay  += "CDR1 AA           :" + CDR1AA + "\n";	
		CDRSdisplay  += "unmutated CDR1 AA :" + CodonTrans.findCodons(unMutatedCDR1, 0) + "\n"+"\n";	
		
		CDRSdisplay += "CDR2 Seq          :" + CDR2+"\n";
		CDR2AA= CodonTrans.findCodons(CDR2, 0);
		CDRSdisplay += "CDR2 AA           :" + CDR2AA + "\n";				
		CDRSdisplay += "unmutated CDR2 AA :" + CodonTrans.findCodons(unMutatedCDR2, 0) + "\n"+"\n";
		
		//Zhilong method
		// CDR3 Amino Acid is wrong, the whole seq# do not have stop codon
		// test with seq substring(2) failed, return different results, means not consider for cut seq,				
		/**
		HashMap<String, Integer> Vdict = getV();
		HashMap<String, Integer> Jdict = getJ(); 
		String myCDR3 = getCDR3seq(IGHVName, IGHJName, IGHVSeq,
					N1Seq, IGHDSeq, N2Seq, IGHJSeq, Vdict, Jdict);
		System.out.println(myCDR3);				
		System.out.println(CodonTrans.findCodons(myCDR3, 0));
		*/
		
		//IMGT find CDR3 by position
		
		CDR3 = getCDR3seq(seqID, VStart);	
		CDRSdisplay +="CDR3 Seq       :"+ CDR3 +"\n";
		CDR3AA= CodonTrans.findCodons(CDR3, 0);
		CDRSdisplay += "CDR3 AA         :" + CDR3AA + "\n";
			
		int myStart = 3-VStart;
		String seqCodon  =  CodonTrans.findCodons(seqID, myStart);
		
		AminoAcid = seqCodon;	
		AminoAcidDisplay = "Whole Protein Sequence: "+ "\n" + AminoAcid +"\n";
		if (! seqCodon.endsWith("-STOP"))  
			AminoAcidDisplay  += "-NO STOP Codon";
		else
			AminoAcidDisplay  += "-STOP Codon";
		
		
		
		} catch (Exception IOE) {
			throw new Error(
					"exception found : " + IOE.getMessage());
		}
		
		
	}
	
	public String getCDRSdisplay(){
		return CDRSdisplay;
	}
	
	public String getCDR1(){
		return CDR1;
	}
	
	
	public String getCDR2(){
		return CDR2;
	}
	
	public String getCDR3(){
		return CDR3;
	}
	
	public String getCDR1AA(){
		return CDR1AA;
	}
	
	public String getCDR2AA(){
		return CDR2AA;
	}
	
	public String getCDR3AA(){
		return CDR3AA;
	}
	
	public String getAminoAcid(){
		return AminoAcid;
	}
	
	public String getAminoAcidDisplay(){
		return AminoAcidDisplay;
	}

	
	private Alignment doAlignment(JAligner jaligner,
			String UMS_string_uppercase, String VGene_string) throws Exception {	
		
		final float OPEN_GAP_COST = (float) 10.0;
		final float EXTEND_GAP_COST = (float) 0.5;
		final String SCORE_MATRIX_NAME = "MATCH";
		
		String VGene_string_uppercase = VGene_string.toUpperCase();
		Alignment alignment_result = JAligner.sw(UMS_string_uppercase,
				VGene_string_uppercase, SCORE_MATRIX_NAME, OPEN_GAP_COST,
				EXTEND_GAP_COST);
		
		return alignment_result;
	}//--doAlignment()
	
	
	
	private String FindVgeneSeq(String VgeneName) {
		
		ArrayList VGenes = new ArrayList(300);
		try {
			FastaReader fr = new FastaReader();
			VGenes = fr.readFile(new File("/srvr/ihmmune/tomcat7/webapps/immuFile/IGHV_Repertoire.fa"));			
		} catch (Exception e) {
			e.printStackTrace();
			throw new Error(e.getMessage());
		}

		
		Sequence temp_sequence;
		String seqString="";

		// for every Sequence in arraylist VGenes
		for (int i = 0; i < VGenes.size(); i++) {
			temp_sequence = (Sequence) VGenes.get(i);
		
			if( temp_sequence.getName().equals(VgeneName)) {	
			seqString = temp_sequence.seqString();
			break;
			}
			
		}//--for(i)

		return seqString;

	}//--FindVgeneSeq()

	
	public String getCDR1seq(String VgeneSeq, int VStart) {
		
		int start= 3*27 - VStart;
		int end= 3*38 - VStart;
		return VgeneSeq.substring(start, end) ;
	}
	
	public String getCDR2seq(String VgeneSeq, int VStart) {
		
		int start= 3*56- VStart;
		int end= 3*65 - VStart;
		return VgeneSeq.substring(start, end) ;
	}
	
//get CDR3seq by position, correct amino acid	
    public String getCDR3seq(String wholeSeq, int VStart) {
		
		int start= 3*105- VStart;
		int end= 3*116 - VStart;
		return wholeSeq.substring(start, end) ;
	}

    //ZhiLong method get CDR3seq, maybe more accurate, 
    //but not identify correct amino acid
    //substring start from 2 return different results, means not consider fragment seq
	public String getCDR3seq(String IGHVName, String IGHJName, String IGHVSeq,
			String N1Seq, String IGHDSeq, String N2Seq, String IGHJSeq,
			HashMap<String, Integer> Vdict, HashMap<String, Integer> Jdict) {
		try {
			
			String vname = IGHVName.split("[ _/(]")[0].replace("\"", "");
		
			if (IGHVName.contains("_")) {
				if (IGHVName.indexOf("_") < IGHVName.indexOf("IGHV")) {
					vname = IGHVName.split("[_/(]")[1].replace("\"", "");
				}
			}
			String jname = "";
			if(IGHJName.contains("[")){
			jname = IGHJName.split("[ /(]")[0].replace("\"", "");
			}
			else{
				jname = IGHJName;
			}
			
			String start = IGHVSeq.replace("\"", "").substring(
					(IGHVSeq.length() - Vdict.get(vname))).replace(".", "");
			String end = IGHJSeq.replace("\"", "").substring(0,
					Jdict.get(jname) + 1).replace(".", "");
			String cdr3seq = start + N1Seq + IGHDSeq.replace(".", "") + N2Seq
					+ end;
		
			return cdr3seq;
		} catch (Exception exception) {
			throw new Error("IGHV not in repertoire:" + IGHVName);
		}
	}
	
   //only in Zhiling methong
	public HashMap<String, Integer> getV() throws Exception {
		
		BufferedReader vreader = new BufferedReader(new InputStreamReader(
				new FileInputStream(new File("src/IGHV_Repertoire_CDR3.csv"))));
		String vdictinput = null;
		HashMap<String, Integer> Vdict = new HashMap<String, Integer>();

		while ((vdictinput = vreader.readLine()) != null) {
			String[] VdictList = vdictinput.split(",");
			String[] VGeneName = VdictList[0].split("[_(]");
			String VGene = VGeneName[0];
			Integer VSeqLen = Integer.parseInt(VdictList[2]);
			Vdict.put(VGene, VSeqLen);
		}
		
		vreader.close();
		//System.out.println(Vdict);
		return Vdict;
	}
	
	//only in Zhiling methong
	public HashMap<String, Integer> getJ() throws Exception {
		
		BufferedReader jreader = new BufferedReader(new InputStreamReader(
				new FileInputStream(new File("src/IGHJ_Repertoire_CDR3.csv"))));
		String jdictinput = null;
		HashMap<String, Integer> Jdict = new HashMap<String, Integer>();
		
		while ((jdictinput = jreader.readLine()) != null) {
			String[] JdictList = jdictinput.split(",");
			String JGene = JdictList[0];
			Integer JSeqLen = Integer.parseInt(JdictList[2]);
			Jdict.put(JGene, JSeqLen);
		}
		
		jreader.close();
		//System.out.println(Jdict);
		return Jdict;
	}

}

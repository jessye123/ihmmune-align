package immu;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;

import iHMMuneAlign.AlignmentThread;
import iHMMuneAlign.MutabilityScoreHotspot;
import iHMMuneAlign.ProbabilityHolder;

public class SingleCellTool {
	
	 IHMMuneResult myIHMMue;
	 String AlignmatchResult;
	 

	public void SingleCellAlign(String mySeq){
		
		File V=   new File("/srvr/ihmmune/tomcat7/webapps/immuFile/IGHV_Repertoire.fa");
		File D =  new File("/srvr/ihmmune/tomcat7/webapps/immuFile/IGHD_Repertoire.fa");
		File J=   new File("/srvr/ihmmune/tomcat7/webapps/immuFile/IGHJ_Repertoire.fa");

	   File openFile = new File("/srvr/ihmmune/tomcat7/webapps/immuFile/Current iHMMune-align Probabilities.PH");
	   
	    
		ProbabilityHolder ph ;		
        try {	       	
			BufferedReader br = new BufferedReader(new FileReader(openFile));		  	
		  	String input = "";		  	
		  	String line = br.readLine();
		  	while(line != null)
		  	{
		  		input += (line + " "); // extra spacing because new line is removed
		  		line = br.readLine();
		  		System.out.println(line);
		  	}
		    // create a new Probability Holder from the input (contents) of the file
		    ph = ProbabilityHolder.createFromString(input);   
		    br.close();
		} catch (IOException ioe) {
			throw new Error("loadProbsButton failed:  " + ioe.getMessage());
		}
        
		
		MutabilityScoreHotspot mu= new MutabilityScoreHotspot();	
		
		/*   input is Single sequence  */
		
		 System.out.println("before run");
		
		try{
	      Sequence s=  DNATools.createDNASequence( mySeq,"");
	    
	       AlignmentThread aThread= new AlignmentThread(null,s, V,D,J,ph,mu,(byte) 1,1);
	       aThread.runAlignment();	
	     
	        myIHMMue = aThread.getIHMMuneResult();
	       AlignmatchResult = aThread.getAlignmatchResult();
	     }
	     catch (Exception e){  }
		 System.out.println("after run");
		
	
	}
	
	public IHMMuneResult getIHMMuneResult (){
		return myIHMMue;
	}
	
	public String getAlignmatchResult (){
		return AlignmatchResult;
	}
	
	
}

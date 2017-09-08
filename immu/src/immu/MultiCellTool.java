package immu;

import java.io.*;


import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;

import iHMMuneAlign.AlignmentThread;
import iHMMuneAlign.MutabilityScoreHotspot;
import iHMMuneAlign.ProbabilityHolder;

public class MultiCellTool {
	
	
	 String AlignResultFile;
	 
	//main for testing Only
	public static void main(String[] args) throws Exception {
      
		MultiCellTool myTool= new MultiCellTool();
		myTool.MultiCellAlign("C:\\Users\\test\\webworkspace\\immu\\src\\Andrews Sequences set one.txt",false);
		System.out.println(myTool.getAlignResultFile());
		
	}

	public void MultiCellAlign(String inputFile, Boolean CDRS){
		
String mySeq;  // the Seq align with
		
		File V= new File("/srvr/ihmmune/tomcat7/webapps/immuFile/IGHV_Repertoire.fa");
		File D = new File("/srvr/ihmmune/tomcat7/webapps/immuFile/IGHD_Repertoire.fa");
		File J= new File("/srvr/ihmmune/tomcat7/webapps/immuFile/IGHJ_Repertoire.fa");

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
		  	}
		    // create a new Probability Holder from the input (contents) of the file
		    ph = ProbabilityHolder.createFromString(input);   
		    br.close();
		} catch (IOException ioe) {
			throw new Error("loadProbsButton failed:  " + ioe.getMessage());
		}
        
		
		MutabilityScoreHotspot mu= new MutabilityScoreHotspot();
	
		System.out.println("before run");
		 try {	
			
			 if(CDRS)  AlignResultFile= "IHMMCDRS"+ inputFile;
			 else     AlignResultFile= "IHMM"+ inputFile;
			 //AlignResultFile= "/srvr/ihmmune/tomcat7/webapps/immuFile/"+ AlignResultFile;		
			 AlignResultFile = AlignResultFile.replaceAll("txt", "csv");
			 AlignResultFile=AlignResultFile+".csv";
	         inputFile= "/srvr/ihmmune/tomcat7/webapps/immuFile/"+ inputFile;
	      	File openFile2 = new File(inputFile);
		    PrintWriter writer = new PrintWriter("/srvr/ihmmune/tomcat7/webapps/immuFile/"+AlignResultFile, "UTF-8");
		    
		  writer.print("#seq"+","+"VGene_name");  writer.print(","+"DGene_name");  writer.print(","+"JGene_name");  
		  writer.print(","+"V_seq"); writer.print(","+"N1_seq");writer.print(","+"D_seq");writer.print(","+"N2_seq");writer.print(","+"J_seq");
		  writer.print(","+"V_mut");   writer.print(","+"D_mut");   writer.print(","+"J_mut"); 
		  if(CDRS) {
			  writer.print(","+"CDR1");   writer.print(","+"CDR2");   writer.print(","+"CDR3"); 
			  writer.print(","+"CDR1_AA");   writer.print(","+"CDR2_AA");   writer.print(","+"CDR3_AA"); 
			  writer.print(","+"Seq_AA");
		  }
		       	
				BufferedReader br2 = new BufferedReader(new FileReader(openFile2));		  			  	
			  	String line = br2.readLine();
			  	while(line != null)   
			  	{		
			  		line = br2.readLine();	  		
			  		if (line!=null && !line.startsWith(">")){		  			
			  	    mySeq=line;
			  	    Sequence s=  DNATools.createDNASequence(mySeq,"");
			 		AlignmentThread aThread= new AlignmentThread(writer,s, V,D,J,ph,mu,(byte) 1,1);
			 		aThread.runAlignment();	
			 		 if(CDRS) {
			 			IHMMuneResult  myIHMMue = aThread.getIHMMuneResult(); 
			 			FindCDRs findCDRs= new FindCDRs();
			 			findCDRs.GetCDRs(myIHMMue );
			 			//PrintWrite write CDRs information here
			 			 writer.print(","+findCDRs.getCDR1());   writer.print(","+findCDRs.getCDR2());   writer.print(","+findCDRs.getCDR3()); 
			 			 writer.print(","+findCDRs.getCDR1AA());   writer.print(","+findCDRs.getCDR2AA());   writer.print(","+findCDRs.getCDR3AA()); 			
						 writer.print(","+findCDRs.getAminoAcid());
			 		 } //end if(CDRS)
			 		
			  		}// end if
			  	}//end while  
			    
			    br2.close(); 
			    writer.close();
			} catch (Exception ioe) {
				throw new Error("loadProbsButton failed:  " + ioe.getMessage());
			}
		
		System.out.println("after run");
		
	
	}
	
	
	public String getAlignResultFile(){
		return AlignResultFile;
	}
	
	/*
	public void addColumn(String path,String fileName) throws IOException{
	    BufferedReader br=null;
	    BufferedWriter bw=null;
	    final String lineSep=System.getProperty("line.separator");

	    try {
	        File file = new File(path, fileName);
	        File file2 = new File(path, fileName+".1");//so the
	                    //names don't conflict or just use different folders

	        br = new BufferedReader(new InputStreamReader(new FileInputStream(file))) ;
	        bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file2)));
	        String line = null;
	                    int i=0;
	        for ( line = br.readLine(); line != null; line = br.readLine(),i++)
	        {               
                 //get IHMMAlign Results, find data
	            String addedColumn = String.valueOf(data.get(i));
	            bw.write(line+addedColumn+lineSep);
	    }

	    }catch(Exception e){
	        System.out.println(e);
	    }finally  {
	        if(br!=null)
	            br.close();
	        if(bw!=null)
	            bw.close();
	    }

	} */
	
}

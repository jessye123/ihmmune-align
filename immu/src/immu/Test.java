package immu;

import iHMMuneAlign.ProbabilityHolder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class Test {
public String myTest(){
	
	return "hello from my class";
}

public String getFile(){
	String filecontent="myFile";
	
	File openFile = new File("/srvr/z2283813/tomcat7/webapps/immuFile/jess.txt");
	   
      try {	       	
			BufferedReader br = new BufferedReader(new FileReader(openFile));		  	
		  	String input = "";		  	
		  	String line = br.readLine();
		  	filecontent+=line;
		  	while(line != null)
		  	{
		  		input += (line + " "); 
		  		line = br.readLine();
		  		
		  	}
		  	filecontent+=input;
		    br.close();
		} catch (IOException ioe) {
			throw new Error("loadProbsButton failed:  " + ioe.getMessage());
		}
      
	
	return filecontent;
}
}

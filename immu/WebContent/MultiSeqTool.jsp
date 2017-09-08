<%@ include file="header.jsp" %>
<%@page import="java.io.*" %>
<%@page import="immu.MultiCellTool, immu.IHMMuneResult, immu.FindCDRs" %>

<nav class="nav">

<form action="/immu/MultiServlet" method="post" enctype="multipart/form-data">
<h2>Multiple Sequences Tool </h2>
<p>Select a fasta File to upLoad:
 <input type="file" name="file">

<p>Or Enter DNA Sequences in fasta format:
<p><textarea rows= "8" cols= "100" name="DNA"></textarea>
 <p style="padding-left: 5em">
 <input type="radio" name="type" value="IHMM" checked="true"> IHMMuneAlign Output
 <p style="padding-left: 5em">
 <input type="radio" name="type" value="IHMMcdr"> IHMMuneAlign with CDRs and Protein Sequencing<br>
 <p style="padding-left: 5em"><input value=" GO " type="submit">
</form>
</nav>

<article class="article">

<% if (request.getAttribute("message") != null) {
	
	String originalFile="http://ihmmune.srvr.cse.unsw.edu.au/immuFile/"+request.getAttribute("filename"); %>

       <P>  <%= request.getAttribute("message")%>
       <P>Right Click DownLoad <a href="<%=originalFile%>" download>The Original File</a>
     
  <h1>Multiple Sequences Analysis Results:</h1> 
  
    <% 
    if ( (request.getAttribute("type")).equals("IHMM")) {
    MultiCellTool mytool= new MultiCellTool(); 
    String filename= ""+request.getAttribute("filename");
    mytool.MultiCellAlign( filename,false );
    String myResultFile = mytool.getAlignResultFile(); 
    myResultFile= "http://ihmmune.srvr.cse.unsw.edu.au/immuFile/"+myResultFile;
    %>  
    <p>IHMMune Alignment Only: 
    <p>Right Click DownLoad <a href="<%=myResultFile%>" download>The Result file</a>
   <%   }
    
    else {    MultiCellTool mytool= new MultiCellTool(); 
    String filename= ""+request.getAttribute("filename");
    mytool.MultiCellAlign( filename, true );
    String myResultFile = mytool.getAlignResultFile();
    myResultFile= "http://ihmmune.srvr.cse.unsw.edu.au/immuFile/"+myResultFile;
    %>  
     <p>IHMMune Alignment With CDRS and Protein Sequencing: 
     <P>Right Click DownLoad <a href="<%=myResultFile%>" download>The Result File</a>
   <%   
    }
    
  } %>

</article>

<%@ include file="footer.jsp" %>
<%@ page import="java.io.*" %>
<%@page import="immu.SingleCellTool, immu.IHMMuneResult, immu.FindCDRs" %>
<%@ include file="header.jsp" %>


<nav class="nav">
<h2>Single Sequence Tool</h2>
<p>Enter IGH DNA Sequence:
<p><form action="SingleSeqTool.jsp" method="post">
<textarea rows= "4" cols= "100" name="DNA">caggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgcccctcacctgcggtgtctatggtgggtccttcactggtgacttctggacctggatccgccagcccccagggaagggactggagtggattggggaaatctatcaaagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcaatagccacgtccaagaaccaattctccctgaggctgaattctttgaccgccgcggacacggccaaatatttctgtgcgagaggcctctcgaatactgcaggtcgtcggggcccacccgctaaggctatggacgtctggggccaagggaccacggtcaccgtctcctca</textarea>
 
<input value=" GO " type="submit">
</form>
</nav>

<article class="article">
<%
// Grab the variables from the form.
  String DNA = request.getParameter("DNA");
  String type = request.getParameter("type");
 
%>
<%-- Processing data. --%>
<% if( DNA != null ) { %>

<h1>Single Sequence Analysis Results:</h1> 

<% 

SingleCellTool mytool= new SingleCellTool(); 
mytool.SingleCellAlign( DNA );
IHMMuneResult myResult = mytool.getIHMMuneResult (); 
String myAlignResult =mytool.getAlignmatchResult();
if(myResult !=null) {%>

<P> DNA Sequence: <pre><%= myResult.getSeqID() %></pre>
<P> Match Display: <P>
<pre><%= myAlignResult %></pre>
<P><pre> V gene name:      <%= myResult.getIGHVName() %> </pre>    
<p><pre> V gene mutation:  <%= myResult.getVMut() %> </pre>
<P><pre> D gene name:      <%= myResult.getIGHDName() %> </pre>
<P><pre> D gene mutation:  <%= myResult.getDMut() %> </pre>
<P><pre> J gene name:      <%= myResult.getIGHJName() %> </pre>   
<P><pre> D gene mutation:  <%= myResult.getDMut() %> </pre>
<% 
//get the CDRS

FindCDRs findCDRs= new FindCDRs();
findCDRs.GetCDRs(myResult);

%>

<P> CDRs informations:
<p><pre><%=findCDRs.getCDRSdisplay() %></pre>

<P> Amino Acid :
<p><pre><%=findCDRs.getAminoAcidDisplay() %></pre>

<%}
}%>

</article>

<%@ include file="footer.jsp" %>
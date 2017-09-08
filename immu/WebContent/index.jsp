<%@ include file="header.jsp" %>


<nav class="nav">
<form action="index.jsp" method="post">
<p style="padding-left: 20em">
 <input type="radio" name="type" value="single" checked="true"> Single Sequence Tool<br>
<p style="padding-left: 20em">
 <input type="radio" name="type" value="multi"> Multiple Sequences Tool<br>
<p style="padding-left: 20em">
 <input type="radio" name="type" value="big"> Graphical Post-analysis Tool<br><br>
<p style="padding-left: 20em">
 <input value=" GO " type="submit">
</form>
</nav>

<article class="article" >
<center>
 <p><img src="antibody.jpg" width=700 height=350>
 
<%
// Grab the variables from the form.  
  String type = request.getParameter("type"); 

%>
<%-- Processing data. --%>
<% if( type != null ) { 
   if(type.equals("single")) { 
	   response.sendRedirect("http://ihmmune.srvr.cse.unsw.edu.au/immu/SingleSeqTool.jsp");
   }  else if(type.equals("multi")){
	   response.sendRedirect("http://ihmmune.srvr.cse.unsw.edu.au/immu/MultiSeqTool.jsp");
   } else {
	   response.sendRedirect("http://ihmmune.srvr.cse.unsw.edu.au/immu/GraphicalTool.jsp");
   }
 } %>
</center>
</article>

<%@ include file="footer.jsp" %>

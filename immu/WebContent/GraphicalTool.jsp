<%@ page import="java.io.*" %>

<%@ include file="header.jsp" %>
<script src="https://ajax.aspnetcdn.com/ajax/jQuery/jquery-3.2.0.min.js"></script>
<script src="http://code.highcharts.com/highcharts.js"></script>
<script src="http://code.highcharts.com/modules/data.js"></script>
<script src="https://code.highcharts.com/highcharts-more.js"></script>
<script src="https://code.highcharts.com/modules/exporting.js"></script>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<script>
function parseCSVFile(allText) {
    var allTextLines = allText.split(/\r\n|\n/);
    var headers = allTextLines[0].split(',');
    var lines = [];

    for (var i=1; i<allTextLines.length; i++) {
        var data = allTextLines[i].split(',');
		for(var j=0;j<data.length;j++){
          if (headers[j]== 'VGene_name' ) {
            lines.push(data[j]);
           }
		}
    }	
	return lines;    
}

function parseCSVVmut(allText) {
    var allTextLines = allText.split(/\r\n|\n/);
    var headers = allTextLines[0].split(',');
    var lines = [];

    for (var i=1; i<allTextLines.length; i++) {
        var data = allTextLines[i].split(',');
		for(var j=0;j<data.length;j++){
        if (headers[j] == 'V_mut') {
            lines.push(data[j]);
        }
		}
    }	
	return lines;    
}

function parseCSVDmut(allText) {
    var allTextLines = allText.split(/\r\n|\n/);
    var headers = allTextLines[0].split(',');
    var lines = [];

    for (var i=1; i<allTextLines.length; i++) {
        var data = allTextLines[i].split(',');
		for(var j=0;j<data.length;j++){
        if (headers[j] == 'D_mut') {
           lines.push(data[j]);
        }
		}
    }	
	return lines;    
}

function parseCSVJmut(allText) {
    var allTextLines = allText.split(/\r\n|\n/);
    var headers = allTextLines[0].split(',');
    var lines = [];

    for (var i=1; i<allTextLines.length; i++) {
        var data = allTextLines[i].split(',');
		for(var j=0;j<data.length;j++){
        if (headers[j] == 'J_mut') {
            lines.push(data[j]);
        }
		}
    }	
	return lines;    
}
</script>

<nav class="nav">
<h2>Graphical Post-analysis Tool</h2>
<p> <I>(NOTE: This tool upload IHMMuneAlign output file and give statistics. </I>
<p> <I>You can also use other csv format files, see below for heading requirements:) </I>
<form action="/immu/BigDataServlet" method="post" enctype="multipart/form-data">
 <p>Select a File to upLoad:
 <input type="file" name="file"><BR>
  <p style="padding-left: 5em">
<input type="checkbox" name="type" value="VgenePie"> V Gene Usage Pie <I> (csv file V gene heading: "VGene_name") </I>
 <p style="padding-left: 5em">
<input type="checkbox" name="type" value="VgeneMutation"> V Gene Mutation Distribution  <I> (csv file V mutation heading: "V_mut") </I>
  <p style="padding-left: 5em">
<input type="checkbox" name="type" value="TotalMutation"> V,D,J Gene Mutation Distribution  <I> (csv file V,D,J mutation heading: "V_mut","D_mut","J_mut") </I>
 <p style="padding-left: 5em">
<input value=" GO " type="submit">
</form>
</nav>

<article class="article">
<% if (request.getAttribute("message") != null) {
	
	String originalFile="http://ihmmune.srvr.cse.unsw.edu.au/immuFile/"+request.getAttribute("filename"); %>

       <P>  <%= request.getAttribute("message")%>
       <P>Right Click DownLoad <a href="<%=originalFile%>">The Original File</a>
     
 
    <%   
   if ( request.getAttribute("VgenePie")!= null) {
 
    %>  
    <p>V gene Usage Pie:  
	<div id="container" style="min-width: 310px; height: 400px; max-width: 600px; margin: 0 auto"></div>
	

<script>
 $.get( "http://ihmmune.srvr.cse.unsw.edu.au/immuFile/outputx9x9x9.csv", function(csvfile) {
 var arr= [];
  arr= parseCSVFile(csvfile);
 var result = { };
   for(i=0;i<arr.length;++i) 
   {
       if(!result[arr[i]])
           result[arr[i]]=0;
       ++result[arr[i]];
   }
   
 var employees = [];
 
 for (var i in result){
  employees.push({name: i,data:result[i]}); 
  
 }
  // Highcharts requires the y option to be set
$.each(employees, function (i, point) {
    point.y = point.data;
});
 
Highcharts.chart('container', {
    chart: {
        plotBackgroundColor: null,
        plotBorderWidth: null,
        plotShadow: false,
        type: 'pie'
    },
    title: {
        text: 'V gene Pie'
    },
    tooltip: {
        pointFormat: '{series.name}: <b>{point.percentage:.1f}%</b>'
    },
    plotOptions: {
        pie: {
            allowPointSelect: true,
            cursor: 'pointer',
            dataLabels: {
                enabled: true,
                format: '<b>{point.name}</b>: {point.percentage:.1f} %',
                style: {
                    color: (Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black'
                }
            }
        }
    },
    series: [{
        name: 'Brands',
        colorByPoint: true,
        data: employees
	    
    }]
});

})
   .fail(function() {
    alert( "CSV error" );
  });
</script>
    
 <%   }
   
   if ( request.getAttribute("VgeneMutation")!= null) {
  
    %>  
    <p>V gene Mutation: 
	
	<div id="containerV" style="min-width: 310px; height: 400px; max-width: 600px; margin: 0 auto"></div>
		
<script>
 $.get( "http://ihmmune.srvr.cse.unsw.edu.au/immuFile/outputx9x9x9.csv", function(csvfile) {
 var v1= [];
  v1= parseCSVVmut(csvfile);
 
 var data = [
  {
    y: v1,
    boxpoints: 'all',
    jitter: 0.3,
    pointpos: -1.8,
	name: 'V Gene Mutation',
    type: 'box'
  }
];

Plotly.newPlot('containerV', data);

})
   .fail(function() {
    alert( "CSV error" );
  });
</script>
   
 <%   }
   
   if ( request.getAttribute("TotalMutation")!= null) {
	
	    %>  
	    <p>Total Mutation: 
		<div id="containerVDJ" style="min-width: 310px; height: 400px; max-width: 600px; margin: 0 auto"></div>
	
	    <script>
 $.get( "http://ihmmune.srvr.cse.unsw.edu.au/immuFile/outputx9x9x9.csv", function(csvfile) {
 var v= [];
  v= parseCSVVmut(csvfile);
 var d= [];
  d= parseCSVDmut(csvfile);
  var j = [];
  j= parseCSVJmut(csvfile);
 var vdj=[];
  for (i=0; i<v.length; i++){
  vdj[i]=parseInt(v[i])+parseInt(d[i])+parseInt(j[i]);
  }
 
  
var trace1 = {
  y: v,
  name: "v gene",
  type: 'box'
};

var trace2 = {
  y: d,
  name: "d gene",
  type: 'box'
};  
var trace3 = {
  y: j,
  name: "j gene",
  type: 'box'
};

var trace4 = {
  y: vdj,
  name: "total",
  type: 'box'
};

var data = [trace1, trace2,trace3,trace4];

Plotly.newPlot('containerVDJ', data);

})
   .fail(function() {
    alert( "CSV error" );
  });


</script>

 <%   }  
  } %>

</article>

<%@ include file="footer.jsp" %>
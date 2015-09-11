<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template name="tomtom-to-html.css">
    <![CDATA[
    /* START INCLUDED FILE "tomtom-to-html.css" */
span.A {
  font-size:20px; color:red
}
span.C {
  font-size:20px; color:blue
}
span.G {
  font-size:20px; color:orange
}
span.T {
  font-size:20px; color:green
}
table.targets td {
  padding: 0 10px;
}
table.preview td {
  padding: 0 10px;
}
table.preview tbody td {
  padding-bottom: 10px;
}
.ac {
  text-align: center;
}
.downloadTd {
  padding-left:20px;
}
div.logo_container {
  position:relative;
  width:99%; 
  height:285px; 
  padding:0px;
  margin:0px;
}
img.logo {
  position:absolute; 
  z-index:2;
  max-width:100%;
  max-height:285px;
}
tr.tspace th, tr.tspace td {
  padding-top: 20px;
}
/* motif list link, first style */
a.ml1 {
 background-color: #FFF; 
}
a.ml2 {
  background-color: #FFF;
}
td.ml {
  line-height: 1.8em;
}

table.match_summary {
  width: 100%;
}

table.match_summary th {
  min-width: 10em;
}

table.match_summary div.help {
  opacity: 0.0;
    -moz-transition: opacity 0.5s;
    -webkit-transition: opacity 0.5s;
    -o-transition: opacity 0.5s;
    transition: opacity 0.5s;
}

table.match_summary tr:hover div.help {
  opacity: 1.0;
}
    /* END INCLUDED FILE "tomtom-to-html.css" */
    ]]>
  </xsl:template>
</xsl:stylesheet>


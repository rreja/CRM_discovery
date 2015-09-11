<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template name="mast-to-html.css">
    <![CDATA[
    /* START INCLUDED FILE "mast-to-html.css" */
.more_arrow {
  font-family:Arial,Serif;
  font-size: larger;
  font-weight: bold;
  text-decoration:none; 
}
div.infobox {
  background-color:#ddddff;
  margin-top: 1.6em;
  margin-bottom: 1em;
}
.sequence {
  font-family:Monospace;
}
.sequence_labels {
  position: relative; 
  width: 100%; 
  height: 2.5em; 
  padding: 0px; 
  margin: 0px;
}
.sequence_label {
  font-family:Serif; 
  position:absolute; 
  z-index:2; 
  height:1em; 
  text-align:center; 
  vertical-align:middle;
}
.sequence_label_top {
  top:0px;
}
.sequence_label_bottom {
  top:1.25em;
}
.inlineTitle {
  display:inline;
}
.block_needle {
  position:absolute;
  z-index:4;
  height:30px; 
  width:1px; 
  top:-2px; 
  background-color:gray;
}
.block_handle {
  position:absolute; 
  z-index:5; 
  height:1.1em; 
  width:3em; 
  top:30px; 
  left:-1.5em; 
  text-align:center; 
  vertical-align:middle;
  background-color: LightGrey; 
  border:3px outset grey;
}
.label_container {
  position:relative;
  width:100%;
  height:25px;
  padding:0px; 
  margin: 0px;
}
table.padded-table td { 
  padding:0px 10px; 
}
table.padded-table th { 
  padding:0px 5px; 
}
td.tnum {
  text-align:right;
}
tr.highlight {
  background:#aaffaa;
}
td.dim {
  color: gray;
}
span.dim {
  color: gray;
}

    /* END INCLUDED FILE "mast-to-html.css" */
    ]]>
  </xsl:template>
</xsl:stylesheet>


<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template name="spamo-to-html.css">
    <![CDATA[
    /* START INCLUDED FILE "spamo-to-html.css" */
      /* ++++++++++++++++++++++++++++++ start spamo.css ++++++++++++++++++++++++++++++++++ */
      div.summary_scroll {width: 99%; height: 300px; overflow: auto; margin: 0 auto; margin:1px;}
      /* specifies positions for the titles of the same and opposite strands for both the graph and spacings */
      div.st_box {position:relative; width:15px; height: 260px;}
      div.st_same {position:absolute; width:130px; height:15px; top:57.5px; left:-60px;}
      div.st_oppo {position:absolute; width:130px; height:15px; top:187.5px; left:-60px;}
      /* the 4 quadrants of significant spacings */
/*
      div.sp_box {position:relative; width: 260px; height: 260px;}
      div.sp_same_up {position:absolute; top:0px; height: 129px; left:0px; width:129px; border: 1px solid black;}
      div.sp_same_down {position:absolute; top:0px; height: 129px; left:130px; width:129px; border-top: 1px solid black; border-right: 1px solid black;}
      div.sp_oppo_up {position:absolute; top:130px; height:129px; left:0px; width:129px; border-left: 1px solid black; border-bottom: 1px solid black;}
      div.sp_oppo_down {position:absolute; top:130px; height: 128px; left:130px; width:128px; border: 1px solid black;}
      /* sp_scroll and sp_fixed are for styling the spacing table.
       * they make use of css 3 to style alternating rows with different colours.
*/
      div.sp_box {position:relative; width: 280px; height: 280px;}
      div.sp_same_up {position:absolute; top:0px; height: 139px; left:0px; width:139px; border: 1px solid black;}
      div.sp_same_down {position:absolute; top:0px; height: 139px; left:140px; width:139px; border-top: 1px solid black; border-right: 1px solid black;}
      div.sp_oppo_up {position:absolute; top:140px; height:139px; left:0px; width:139px; border-left: 1px solid black; border-bottom: 1px solid black;}
      div.sp_oppo_down {position:absolute; top:140px; height: 138px; left:140px; width:138px; border: 1px solid black;}
      /* sp_scroll and sp_fixed are for styling the spacing table.
       * they make use of css 3 to style alternating rows with different colours.
       * sp_scroll is used on tables with more than 5 elements to create vertical scrollbars 
       * on the body of the spacing table
       * sp_fixed is for tables with 5 or less elements. It exists because firefox 3 likes
       * to spread the entries out when there are less than needed to fit in the vertical 
       * space*/
      div.sp_scroll {width: 99%; height: 130px; overflow: auto; margin: 0 auto; margin:1px;}
      div.sp_fixed {width: 99%; margin: 0 auto; margin:1px;}
      div.sp_scroll table, div.sp_fixed table {width: 100%; border: none; border-spacing: 0px; background-color: transparent;}
      /* make it scroll for a table body height of larger than 105px, disable the horizontal scrollbar */
      div.sp_scroll table>tbody {overflow: auto; height: 105px; overflow-x: hidden;}
      /* detach table header by making relative and forcing to the top */
      div.sp_scroll table>thead tr {position:relative; top: 0px;}
      /* set alternating row colours (note uses css3 feature nth-child so may not work many places) */
      div.sp_scroll table>tbody tr, div.sp_fixed table>tbody tr {height: auto; background: #ccccee;}
      div.sp_scroll table>tbody tr:nth-child(odd), div.sp_fixed table>tbody tr:nth-child(odd) {background: #eeeeee;}
      /* check box scrolling */
      div.ck_scroll {width: 99%; height: 120px; overflow: auto; margin: 0 auto; margin:1px;}
      div.ck_fixed {width: 99%; margin: 0 auto; margin:1px;}
      /* set column alignment */
      td.sp_gap, td.sp_count {text-align:right;}
      td.sp_gap {padding-right: 0.5em;}
      td.sp_pvalue {text-align:left;}
      table.details>tbody th {text-align:right;}
      /* consensus sequence colours */
      div.seq {font-family: monospace; font-size:20px;}
      div.seqbox {height:61px; margin-right:10px; border: 1px solid black; display: table-cell; vertical-align: middle; text-align: center; padding: 0px 5px;}
      span.A {color:red;}
      span.C {color:blue;}
      span.G {color:orange;}
      span.T {color:green;}
      span.trim {background-color:gray;}
      .deg90rotate {
     -moz-transform: rotate(-90deg);  /* FF3.5+ */
       -o-transform: rotate(-90deg);  /* Opera 10.5 */
  -webkit-transform: rotate(-90deg);  /* Saf3.1+, Chrome */
          transform: rotate(-90deg);  
             filter:  progid:DXImageTransform.Microsoft.Matrix(sizingMethod='auto expand', /* IE6,IE7 */ 
                      M11=6.123233995736766e-17, M12=1, M21=-1, M22=6.123233995736766e-17); 
         -ms-filter: "progid:DXImageTransform.Microsoft.Matrix(M11=6.123233995736766e-17, M12=1, M21=-1, M22=6.123233995736766e-17, sizingMethod='auto expand')"; /* IE8 */
               zoom: 1;

      }
      .sm_group{
        margin: 0px;
        padding: 0px;
        border: 0px;
      }
      .sm_group p{
        font-weight: bold;
        padding-bottom: 0px;
        margin-bottom: 0px;
      }
      .sm_group label{
        display: block;
      }

      .details h4 {
        border: 0px;
      }

      div.pop_content {
        position:absolute;
        z-index:1;
        width:300px;
        padding: 5px;
        background: #E4ECEC;
        font-size: 12px;
        font-family: Arial;
        border-style: double;
        border-width: 3px;
        border-color: #AA2244;
        display:none;
      }

      table.preview {
        border-spacing: 10px 0px;
      }
      
      table.preview tr th {
        text-align: left;
        font-size: 1.2em; 
        font-style: normal;
        font-family: Georgia, "Times New Roman", Times, serif; 
        border-bottom: 1px solid #CCCCCC; 
      }

      table.preview tr td.summarylist {
        line-height: 1.8em;
      }

      /* ++++++++++++++++++++++++++++++ end spamo.css ++++++++++++++++++++++++++++++++++++ */
    /* END INCLUDED FILE "spamo-to-html.css" */
    ]]>
  </xsl:template>
</xsl:stylesheet>


<xsl:stylesheet version='1.0'
	type="text/xsl" 
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns="http://www.w3.org/1999/xhtml"
	xmlns:api="http://bioinformatics.pzr.uni-rostock.de/~guest01/apibaso"
>
 <xsl:output method="html" encoding="ISO-8859-1"/>
 <xsl:template match="/" >
  <html>
   <head>
    <title>Apibaso: apical and basolateral proteins in MDCK cell line</title>
   </head>
   <body>
    <h1 align="center">Apical and Basolateral Sorted Proteins</h1>
     <xsl:apply-templates select="api:variantlist" />
    <hr/>
   </body>
  </html>
 </xsl:template>

 <xsl:template match="api:variantlist">
    <table width="80%" border="1" cellpadding="2" align="center" >
	    <tr bgcolor="orange"><th>Gene / Gene Product (tissue)</th><th>References</th><th>Location</th><th>responsable sequence</th></tr>
    <xsl:apply-templates select="api:variant"/>
    </table>
 </xsl:template>
 	<!-- <xsl:if test="not($uniprot='none')&$empty_string=$uniprot">  -->

 <xsl:template match="api:variant">
   <xsl:variable name="uniprot" select="@UniProt" />
   <xsl:variable name="evidence" select="@evidence" />
   <xsl:variable name="tissue" select="@tissue" />
  <tr>
      <td bgcolor="yellow" valign="top" align="center" rowspan="3">
      <xsl:value-of select="@gene" />
	<xsl:variable name="empty_string" select="''" />
	<br/>
 	<xsl:if test="$empty_string!=$uniprot"> 
		<small><a href="http://www.expasy.org/cgi-bin/niceprot.pl?{$uniprot}">UniProt</a></small>
	</xsl:if> 
 	<xsl:if test="$empty_string!=$evidence"> 
		<small><p>Evidence: <xsl:value-of select="@evidence" /></p></small>
	</xsl:if> 
	<br />

	<xsl:if test="$empty_string!=$tissue"> 
		<small><p>tissue: <xsl:value-of select="@tissue" /></p></small>
	</xsl:if> 
	<br />

      </td>
    
      <td align="center" rowspan="2" >
	<table cellpadding="5">
          <xsl:apply-templates select="api:reference"/>
	</table>
      </td>
      <td valign="top" align="center"><xsl:value-of select="@location" /></td> 
      <td valign="top" align="center"><xsl:apply-templates select="api:sequence-feature"/></td> 
  </tr>
  <tr>
      <td valign="top" align="left" colspan="3"><xsl:value-of select="@synonyms" /></td>
  </tr>
  <tr>
	  <td valign="top" align="left" colspan="3"><small><i><xsl:value-of select="@comment" /></i></small></td> 
  </tr>
 </xsl:template>

 <xsl:template match="api:reference">
  <xsl:variable name="url" select="@href" />
  <tr>
   <td>
      <xsl:variable name="type" select="@type" />
      <xsl:choose> 
       <xsl:when test="'pubmed'=$type"> 
        <xsl:variable name="pmid" select="@pmid" />
    	<a href='http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=pubmed&amp;dopt=Abstract&amp;list_uids={$pmid}'>
        <xsl:value-of select="@title" /></a>
       </xsl:when>
       <xsl:otherwise>
         <a href="{$url}">
         <xsl:value-of select="@title" /></a>
       </xsl:otherwise>
      </xsl:choose>
   </td>
  </tr>
 </xsl:template>

 <xsl:template match="api:sequence-feature">
	<xsl:variable name="sequence" select="@sequence" />
	<xsl:variable name="empty_string" select="' '" />
	
	<xsl:choose>
	<xsl:when test="$empty_string!=$sequence">
		<small><p>found a feature: <xsl:value-of select="@sequence" /></p></small>
	</xsl:when>
	<xsl:otherwise>
         <small><p>found no special sequence! </p></small>
       </xsl:otherwise>
	</xsl:choose>
</xsl:template>
</xsl:stylesheet>

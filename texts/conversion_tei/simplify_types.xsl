<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:math="http://www.w3.org/2005/xpath-functions/math"
    xmlns:tei="http://www.tei-c.org/ns/1.0"
    exclude-result-prefixes="xs math"
    version="3.0">
    
    <!-- To check apps that should be afterwards merged (in XPath 3.1)
    //app[
  let $curr := .
  return 
    following-sibling::node()
      [not(self::comment() or self::processing-instruction())]
      [1]
      [self::app or (self::text() and normalize-space(.) = '')]
      /self::text()
      /following-sibling::node()
        [not(self::comment() or self::processing-instruction())]
        [1]
        [self::app]
        [@type = $curr/@type]
        [count(rdg) = count($curr/rdg)]
        [
          let $next := .
          return every $r1 in $curr/rdg satisfies 
            some $r2 in $next/rdg satisfies 
              deep-equal(
                sort(tokenize($r1/@wit, '\s+')), 
                sort(tokenize($r2/@wit, '\s+'))
              )
        ]
]
    
    -->
    
    
    <xsl:template match="node() | @*">
        <xsl:copy>
            <xsl:apply-templates select="node() | @*"/>
        </xsl:copy>
    </xsl:template>
    
    
    <xsl:template match="tei:app">
        <xsl:copy>
            <xsl:if test="count(tei:rdg) > 1">
                <xsl:attribute name="type">
                <xsl:choose>
                    <xsl:when test="not(@type) and tei:rdg[ normalize-space(.) = '' ]">
                        <xsl:text>sequence</xsl:text>
                    </xsl:when>
                    <xsl:when test="@type = ('graph')">
                        <xsl:text>graphemic</xsl:text>
                    </xsl:when>
                    <xsl:when test="(@type = ('flex', 'morsynt')) or (@type = 'semlex:minor' and @subtype = 'gramm')">
                        <xsl:text>grammatical</xsl:text>
                    </xsl:when>
                    <xsl:when test="(@type = 'semlex:minor' and not(@subtype = 'gramm')) or @type=('semlex:major')">
                        <xsl:text>semantic</xsl:text>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:text>ERROR</xsl:text>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:attribute>
            </xsl:if>
            <xsl:apply-templates select="node() | @*[not(local-name(.) = 'type')]"/>
        </xsl:copy>
    </xsl:template>
    
</xsl:stylesheet>
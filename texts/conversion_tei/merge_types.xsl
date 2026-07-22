<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns="http://www.tei-c.org/ns/1.0"
    xpath-default-namespace="http://www.tei-c.org/ns/1.0"
    xmlns:fn="local-functions"
    exclude-result-prefixes="xs fn">
    
    <xsl:output method="xml" indent="yes"/>
    <xsl:mode on-no-match="shallow-copy"/>
    
    <!-- 2. Helper function to generate a signature for each app.
            Apps with identical @type, rdg count, and sorted @wit sets 
            will produce the exact same signature string. -->
    <xsl:function name="fn:app-signature" as="xs:string">
        <xsl:param name="app" as="element(app)"/>
        
        <xsl:variable name="type" select="string($app/@type)"/>
        
        <!-- Build a canonical signature for each rdg by normalizing and sorting witness sigla -->
        <xsl:variable name="rdg-signatures" as="xs:string*">
            <xsl:for-each select="$app/rdg">
                <xsl:value-of select="string-join(sort(tokenize(normalize-space(@wit), '\s+')), ' ')"/>
            </xsl:for-each>
        </xsl:variable>
        
        <!-- CRITICAL FIX: Sort $rdg-signatures so child <rdg> order doesn't affect the key -->
        <xsl:variable name="sorted-rdg-sigs" select="sort($rdg-signatures)"/>
        
        <!-- Combine type, rdg count, and sorted rdg signatures into one canonical key -->
        <xsl:value-of select="string-join(($type, count($app/rdg), $sorted-rdg-sigs), '||')"/>
    </xsl:function>
    
    <!-- 3. Group child elements inside any parent containing <app> -->
    <xsl:template match="*[app]">
        <xsl:copy>
            <xsl:apply-templates select="@*"/>
            
            <!-- Group adjacent nodes. If a text node between apps is only whitespace, 
                 it receives the adjacent app's group key so it doesn't break grouping. -->
            <xsl:for-each-group select="node()" group-adjacent="
                if (self::app) then 
                fn:app-signature(.)
                else if (self::text() and normalize-space(.) = '' 
                and preceding-sibling::node()[1]/self::app 
                and following-sibling::node()[1]/self::app) then
                fn:app-signature(following-sibling::app[1])
                else 
                generate-id(.)">
                
                <xsl:choose>
                    <!-- If the group contains 2 or more <app> elements, merge them -->
                    <xsl:when test="count(current-group()[self::app]) > 1">
                        <!-- Log the merged/disappeared IDs to the console -->
                        <xsl:message>
                            <xsl:text>Merged group starting with </xsl:text>
                            <xsl:value-of select="current-group()[self::app][1]/@xml:id"/>
                            <xsl:text> | Dropping IDs: </xsl:text>
                            <xsl:value-of select="current-group()[self::app][position() > 1]/@xml:id" separator=", "/>
                        </xsl:message>
                        <xsl:call-template name="merge-apps">
                            <xsl:with-param name="apps" select="current-group()[self::app]"/>
                        </xsl:call-template>
                    </xsl:when>
                    
                    <!-- Otherwise, process nodes normally -->
                    <xsl:otherwise>
                        <xsl:apply-templates select="current-group()"/>
                    </xsl:otherwise>
                </xsl:choose>
                
            </xsl:for-each-group>
        </xsl:copy>
    </xsl:template>
    
    <!-- 4. Logic to combine readings from merged <app>s -->
    <xsl:template name="merge-apps">
        <xsl:param name="apps" as="element(app)*"/>
        <xsl:variable name="first-app" select="$apps[1]"/>
        
        <app type="{$first-app/@type}">
            <xsl:copy-of select="$first-app/@xml:id"/>
            
            <!-- Reconstruct readings for all unique witness configurations across the group -->
            <xsl:for-each select="$first-app/rdg">
                <xsl:variable name="curr-wit-sorted" select="sort(tokenize(normalize-space(@wit), '\s+'))"/>
                
                <rdg wit="{@wit}">
                    <!-- Preserve any extra attributes (like xml:id or corresp) from the first rdg -->
                    <xsl:copy-of select="@*[name() != 'wit']"/>
                    
                    <xsl:for-each select="$apps">
                        <xsl:variable name="matching-rdg" select="rdg[deep-equal(sort(tokenize(normalize-space(@wit), '\s+')), $curr-wit-sorted)]"/>
                        <xsl:apply-templates select="$matching-rdg/node()"/>
                        
                        <xsl:if test="position() != last() and $matching-rdg/node()">
                            <xsl:text> </xsl:text>
                        </xsl:if>
                    </xsl:for-each>
                </rdg>
            </xsl:for-each>
        </app>
    </xsl:template>
    
</xsl:stylesheet>
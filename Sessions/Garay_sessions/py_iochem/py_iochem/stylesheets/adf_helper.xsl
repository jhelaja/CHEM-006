<?xml version="1.0" encoding="UTF-8"?>
<!--

    Create module - Create module inside the ioChem-BD software.
    Copyright © 2014 ioChem-BD (contact@iochem-bd.org)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

-->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:adf="http://www.scm.com/ADF/"
    xmlns:cml="http://www.xml-cml.org/schema"
    xmlns:cmlx="http://www.xml-cml.org/schema/cmlx"
    xmlns:helper="http://www.w3.org/1999/XSL/Helper-Functions"
    exclude-result-prefixes="xs"
    version="2.0">
    
    <!-- Calculation type related variables -->
    <xsl:variable name="adf:GeometryOptimization" select="'Geometry optimization'" />
    <xsl:variable name="adf:SinglePoint" select="'Single point'" />
    <xsl:variable name="adf:TransitionState" select="'Transition state'" />
    <xsl:variable name="adf:Frequencies" select="'Frequencies'" />    
    <xsl:variable name="adf:Quild" select="'Quild'" />    
    <xsl:variable name="adf:NMR" select="'NMR'" />
    <xsl:variable name="adf:GeometryConverged" select="'(Geometry Converged|GEOMETRY CONVERGED|CONVERGED)'"/>
    <xsl:variable name="adf:GeometryNotConverged" select="'GEOMETRY DID NOT CONVERGE'"/>        
    
    <xsl:function name="adf:getCalcTypes">
        <xsl:param name="runType" as="node()*"/>
        <xsl:param name="isNMR" as="xs:boolean"/>
        
        <xsl:variable name="calcType">
            <xsl:for-each select="$runType">
                <xsl:choose>
                    <xsl:when test="compare(.,'GEOMETRY OPTIMIZATION') = 0">
                        <xsl:value-of select="$adf:GeometryOptimization"/>                    
                    </xsl:when>
                    <xsl:when test="compare(.,'SINGLE POINT') = 0">
                        <xsl:value-of select="$adf:SinglePoint"/> 
                    </xsl:when>
                    <xsl:when test="compare(.,'TRANSITION STATE') = 0">
                        <xsl:value-of select="$adf:TransitionState"/>
                    </xsl:when>
                    <xsl:when test="compare(.,'FREQUENCIES') = 0">
                        <xsl:value-of select="$adf:Frequencies"/>
                    </xsl:when>                    
                </xsl:choose> 
                <xsl:value-of select="' '"/>
            </xsl:for-each>
        </xsl:variable>
        
        <xsl:variable name="nmr" select="
            if($isNMR)
            then concat(' ',$adf:NMR)
            else
            ''"
        />
        <xsl:sequence select="concat(adf:trim($calcType), $nmr)"/>                 
    </xsl:function>
    
    
    <xsl:function name="adf:getCalcType">
        <xsl:param name="runType" as="xs:string"/>
        <xsl:param name="isQuild" as="xs:boolean"/>
        <xsl:param name="isNMR" as="xs:boolean"/>
        <!-- Need to add LINEARTRANSIT, INTRINSICREACTIONCOORDINATE types  -->
        <xsl:variable name="calcType" select="
            if(compare($runType,'GEOMETRY OPTIMIZATION') = 0) 
            then $adf:GeometryOptimization 
            else if(compare($runType,'SINGLE POINT') = 0)
            then $adf:SinglePoint
            else if(compare($runType,'TRANSITION STATE') = 0)
            then $adf:TransitionState
            else if(compare($runType,'FREQUENCIES') = 0)
            then $adf:Frequencies
            else
            $adf:SinglePoint"
        />                                                         
        <xsl:variable name="quild" select="
            if($isQuild)
            then concat(' ',$adf:Quild)
            else
            ''"
        />
        <xsl:variable name="nmr" select="
            if($isNMR)
            then concat(' ',$adf:NMR)
            else
            ''"
        />
        <xsl:sequence select="concat($calcType, $quild, $nmr)"/>         
    </xsl:function>
    
    
    <!-- Correct intensities on SCANFREQ -->
    <xsl:function name="adf:getIntensities">
        <xsl:param name="intensities" />
        <xsl:param name="scanfreq"  />
        
        <xsl:choose>
            <xsl:when test="not(exists($scanfreq))">
                <xsl:copy-of select="$intensities" />
            </xsl:when>
            <xsl:otherwise>
                <!-- Generate translation map -->
                <xsl:variable name="translations">                                        
                    <xsl:variable name="oldfreqs" select="tokenize($scanfreq/cml:list/cml:array[@dictRef='a:oldfreq'], '\s+')"/>
                    <xsl:variable name="newfreqs" select="tokenize($scanfreq/cml:list/cml:array[@dictRef='a:newfreq'], '\s+')"/>
                    <xsl:variable name="intensities" select="tokenize($scanfreq/cml:list/cml:array[@dictRef='cc:irintensity'], '\s+')"/>                    
                    <xsl:for-each select="1 to count($oldfreqs)">                            
                        <xsl:variable name="outerIndex" select="."/>                                                                                                                                                
                        <xsl:element name="frequency">
                            <xsl:attribute name="oldfreq" select="$oldfreqs[$outerIndex]" />
                            <xsl:attribute name="newfreq" select="$newfreqs[$outerIndex]" />
                            <xsl:attribute name="intensities" select="$intensities[$outerIndex]" />
                        </xsl:element>
                    </xsl:for-each>                    
                </xsl:variable>
                
                <!-- Build the arrays with replaced values -->      
                <xsl:variable name="sourceFrequencies" select="tokenize($intensities/cml:array[@dictRef='cc:frequency'], '\s+')"/>
                <xsl:variable name="finalFrequencies">
                    <xsl:for-each select="1 to count($sourceFrequencies)">
                        <xsl:variable name="outerIndex" select="." />
                        <xsl:variable name="value" select="round-half-to-even(number($sourceFrequencies[$outerIndex]), 3)" />
                        <xsl:value-of select="if(exists($translations/frequency[@oldfreq = $value])) then 
                            concat($translations/frequency[@oldfreq = $value]/@newfreq, ' ')
                            else
                            concat($value, ' ')
                            "/>
                    </xsl:for-each>
                </xsl:variable>                
                
                <xsl:variable name="sourceIntensities" select="tokenize($intensities/cml:array[@dictRef='cc:absortion'], '\s+')"/>
                <xsl:variable name="finalIntensities">
                    <xsl:for-each select="1 to count($sourceIntensities)">
                        <xsl:variable name="outerIndex" select="." />
                        <xsl:variable name="value" select="round-half-to-even(number($sourceFrequencies[$outerIndex]), 3)" />
                        <xsl:value-of select="if(exists($translations/frequency[@oldfreq = $value])) then 
                            concat($translations/frequency[@oldfreq = $value]/@intensities, ' ')
                            else
                            concat($value, ' ')
                            "/>
                    </xsl:for-each>
                </xsl:variable>
                
                <xsl:variable name="sourceDipoles" select="tokenize($intensities/cml:array[@dictRef='cc:dipole'], '\s+')"/>
                <xsl:variable name="finalDipoles">
                    <xsl:for-each select="1 to count($sourceDipoles)">
                        <xsl:variable name="outerIndex" select="." />
                        <xsl:variable name="value" select="round-half-to-even(number($sourceFrequencies[$outerIndex]), 3)" />
                        <!-- Calculate dipole stength as intensity / (freq * 2.506643e-04) -->
                        <xsl:value-of select="if(exists($translations/frequency[@oldfreq = $value])) then 
                            
                            concat(
                            format-number(
                            number($translations/frequency[@oldfreq = $value]/@intensities) 
                            div (number($translations/frequency[@oldfreq = $value]/@newfreq) * 2.506643e-04),
                            '#0.000')
                            , ' ')
                            else
                            concat($value, ' ')
                            "/>
                    </xsl:for-each>
                </xsl:variable>
                
                <!-- Return the corrected module -->
                <xsl:element name="module" namespace="http://www.xml-cml.org/schema">
                    <xsl:attribute name="cmlx:templateRef" select="'intensities'"/>
                    <xsl:element name="array" inherit-namespaces="no" namespace="http://www.xml-cml.org/schema">
                        <xsl:attribute name="dataType" select="'xsd:double'"/>
                        <xsl:attribute name="dictRef" select="'cc:frequency'"/>
                        <xsl:attribute name="units" select="'nonsi:cm-1'" />
                        <xsl:attribute name="size" select="$intensities/cml:array[@dictRef='cc:frequency']/@size"/>                       
                        <xsl:value-of select="adf:trim($finalFrequencies)"/>
                    </xsl:element>
                    
                    <xsl:element name="array" inherit-namespaces="no" namespace="http://www.xml-cml.org/schema">
                        <xsl:attribute name="dataType" select="'xsd:double'"/>
                        <xsl:attribute name="dictRef" select="'cc:dipole'"/>
                        <xsl:attribute name="units" select="'nonsi2:1e-40.esu2.cm2'" />
                        <xsl:attribute name="size" select="$intensities/cml:array[@dictRef='cc:dipole']/@size"/>                       
                        <xsl:value-of select="adf:trim($finalDipoles)"/>
                    </xsl:element>
                    
                    <xsl:element name="array" inherit-namespaces="no" namespace="http://www.xml-cml.org/schema">
                        <xsl:attribute name="dataType" select="'xsd:double'"/>
                        <xsl:attribute name="dictRef" select="'cc:absortion'"/>
                        <xsl:attribute name="units" select="'nonsi:cm-1'" />
                        <xsl:attribute name="size" select="$intensities/cml:array[@dictRef='cc:absortion']/@size"/>                       
                        <xsl:value-of select="adf:trim($finalIntensities)"/>
                    </xsl:element>
                </xsl:element>                   
                
            </xsl:otherwise>
        </xsl:choose>       
    </xsl:function>
    
    <!-- Correct frequencies on SCANFREQ -->
    <xsl:function name="adf:getFrequencies">
        <xsl:param name="frequencies" />
        <xsl:param name="scanfreq"  />
        
        <xsl:choose>
            <xsl:when test="not(exists($scanfreq))">
                <xsl:copy-of select="$frequencies" />
            </xsl:when>
            <xsl:otherwise>
                <!-- Generate translation map -->
                <xsl:variable name="translations">                                        
                    <xsl:variable name="oldfreqs" select="tokenize($scanfreq/cml:list/cml:array[@dictRef='a:oldfreq'], '\s+')"/>
                    <xsl:variable name="newfreqs" select="tokenize($scanfreq/cml:list/cml:array[@dictRef='a:newfreq'], '\s+')"/>                
                    <xsl:for-each select="1 to count($oldfreqs)">                            
                        <xsl:variable name="outerIndex" select="."/>                                                                                                                                                
                        <xsl:element name="frequency">
                            <xsl:attribute name="oldfreq" select="$oldfreqs[$outerIndex]" />
                            <xsl:attribute name="newfreq" select="$newfreqs[$outerIndex]" />                       
                        </xsl:element>
                    </xsl:for-each>                    
                </xsl:variable>
                <!-- Build an array with replaced values -->                
                <xsl:variable name="sourceFrequencies" select="tokenize($frequencies/cml:array[@dictRef='cc:frequency'], '\s+')"/>
                <xsl:variable name="finalFrequencies">
                    <xsl:for-each select="1 to count($sourceFrequencies)">
                        <xsl:variable name="outerIndex" select="." />
                        <xsl:variable name="value" select="round-half-to-even(number($sourceFrequencies[$outerIndex]), 3)" />
                        <xsl:value-of select="if(exists($translations/frequency[@oldfreq = $value])) then 
                            concat($translations/frequency[@oldfreq = $value]/@newfreq, ' ')
                            else
                            concat($value, ' ')
                            "/>
                    </xsl:for-each>
                </xsl:variable>
                
                <!-- Return the corrected module -->
                <xsl:element name="module" namespace="http://www.xml-cml.org/schema">
                    <xsl:attribute name="cmlx:templateRef" select="'vibrations'"/>
                    <xsl:attribute name="dictRef" select="'cc:vibrations'"/>
                    <xsl:element name="array" inherit-namespaces="no" namespace="http://www.xml-cml.org/schema">
                        <xsl:attribute name="dataType" select="'xsd:double'"/>
                        <xsl:attribute name="dictRef" select="'cc:frequency'"/>
                        <xsl:attribute name="units" select="'nonsi:cm-1'" />
                        <xsl:attribute name="size" select="$frequencies/cml:array[@dictRef='cc:frequency']/@size"/>                       
                        <xsl:value-of select="adf:trim($finalFrequencies)"/>
                    </xsl:element>
                    <xsl:copy-of select="$frequencies/cml:array[@dictRef='cc:elementType']" />
                    <xsl:copy-of select="$frequencies/cml:array[@dictRef='cc:displacement']" />
                </xsl:element>                                  
            </xsl:otherwise>
        </xsl:choose>       
    </xsl:function>
    
    
    
    <!-- Trim whitespaces from string -->
    <xsl:function name="adf:trim" as="xs:string">
        <xsl:param name="arg" as="xs:string?"/>        
        <xsl:sequence select="replace(replace($arg,'\s+$',''),'^\s+','')"/>        
    </xsl:function>
    
</xsl:stylesheet>
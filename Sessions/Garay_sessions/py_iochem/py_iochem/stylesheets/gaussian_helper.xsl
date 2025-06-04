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
    xmlns:gaussian="http://www.gaussian.com/"
    exclude-result-prefixes="xs"
    version="2.0">
    
    
    <!-- Calculation type related variables -->
    <xsl:variable name="gaussian:GeometryOptimization" select="'Geometry optimization'" />
    
    <xsl:variable name="shellWaveFunctions" select="'RO|U|R|'"/>
    <xsl:variable name="shellWaveFunctionsObligatory" select="'RO|U|R'"/>
    <xsl:variable name="methodsRegex" select="upper-case('DFT|MM|Amber|Dreiding|UFF|AM1|PM3|PM3MM|PM6|PDDG|HF|HFS|XAlpha|HFB|VSXC|HCTH|HCTH93|HCTH147|HCTH407|tHCTH|M06L|B97D|LSDA|LC-wPBE|CAM-B3LYP|wB97XD|wB97|wB97X|MN15|M11|SOGGA11X|N12SX|MN12SX|PW6B95|PW6B95D3|M08HX|LC-BLYP|B3LYP|B3P86|B3PW91|B1B95|mPW1PW91|mPW1LYP|mPW1PBE|mPW3PBE|B98|B971|B972|PBE1PBE|B1LYP|O3LYP|TPSSh|BMK|M06|M06HF|M062X|M05|M052X|HISSbPBE|X3LYP|BHandH|BHandHLYP|tHCTHhyb|HSEh1PBE|HSE2PBE|HSEhPBE|OHSE2PBE|OHSE1PBE|PBEh1PBE|CASSCF|CAS|MP2|MP3|MP4|MP5|B2PLYP|B2PLYPD|mPW2PLYP|mPW2PLYPD|QCISD|CCD|CCSD\(T\)|CCSD|CC|QCID|BD|EPT|CBS-4M|CBS-QB3|ROCBS-QB3|CBS-APNO|G1|G2|G2MP2|G3|G3MP2|G3B3|G3MP2B3|G4|G4MP2|W1U|W1BD|W1RO|CIS|CIS\(D\)|CID|CISD|TD|EOMCCSD|ZINDO|DFTB|DFTBA|GVB|CNDO|INDO|MINDO|MNDO|NMR|SAC-CI')"/>
    <xsl:variable name="multipleMethodRegex" select="'ONIOM|IRCMAX'"/>
    <xsl:variable name="exchangeFunctional" select="upper-case('S|XA|B|PW91|mPW|G96|PBE|OPBE|O|TPSS|BRx|PKZB|wPBEh|PBEh|LC-')"/>
    <xsl:variable name="correlationFunctional" select="upper-case('VWN|VWN5|LYP|PL|P86|PW91|B95|PBE|TPSS|KCIS|BRC|PKZB|VP86|V5LYP')"/>
    
    
    <xsl:function name="gaussian:getCalcType" as="xs:string">        
        <xsl:param name="isOptimization" as="xs:boolean"/>
        <xsl:param name="hasStationaryPoint" as="xs:boolean"/>
        <xsl:param name="hasMinimum" as="xs:boolean"/>        
        <xsl:value-of select="gaussian:getCalcType($isOptimization, $hasStationaryPoint, $hasMinimum, false())"/>        
    </xsl:function>
    
    <xsl:function name="gaussian:getCalcType" as="xs:string">        
        <xsl:param name="isOptimization" as="xs:boolean"/>
        <xsl:param name="hasStationaryPoint" as="xs:boolean"/>
        <xsl:param name="hasMinimum" as="xs:boolean"/>
        <xsl:param name="isEET" as="xs:boolean"/>
        
        <xsl:choose>
            <xsl:when test="$isEET">
                <xsl:sequence select="'Excitation energy transfer'"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:variable name="calcType" select="if($isOptimization) then 'Geometry optimization' else 'Single point'"/>       
                <xsl:choose>
                    <xsl:when test="$hasStationaryPoint">
                        <xsl:variable name="hasMinimum" select="if($hasMinimum) then ' Minimum' else ' TS'"/>
                        <xsl:sequence select="concat($calcType, ' ' , $hasMinimum)"/>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:sequence select="concat($calcType, ' Structure')"/>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:function>
    
    <!-- If we must append new gaussian methods or functionals, we must mind to escape conflictive characters as ( ) . * -->
    
    <xsl:function name="gaussian:isMethod">
        <xsl:param name="candidate" />
        
        
        <xsl:variable name="methodRegex" select="concat('(', $shellWaveFunctions,')(',$methodsRegex ,')')"/>
        <xsl:variable name="compoundMethodRegex" select="concat('(', $shellWaveFunctions, ')(', $exchangeFunctional, ').*(', $correlationFunctional, ')')"/>
        
        <xsl:variable name="isMethod" select="if(matches(upper-case($candidate),$multipleMethodRegex)
            or matches(upper-case($candidate),$methodRegex) 
            or matches(upper-case($candidate),$compoundMethodRegex)) 
            then true()
            else false()"/>
        <xsl:sequence select="$isMethod"/>        
    </xsl:function>
    
    <xsl:function name="gaussian:methodSpinType">
        <xsl:param name="candidate" />
        
        <xsl:variable name="methodRegex" select="concat('(', $shellWaveFunctionsObligatory,')(',$methodsRegex ,')')"/>
        <xsl:variable name="compoundMethodRegex" select="concat('(', $shellWaveFunctionsObligatory, ')(', $exchangeFunctional, ').*(', $correlationFunctional, ')')"/>        
        <xsl:variable name="methodSpinType" select="if(matches(upper-case($candidate),$methodRegex) 
            or matches(upper-case($candidate),$compoundMethodRegex)) 
            then
            if(starts-with($candidate,'R'))                    
            then 'Restricted' 
            else    
            'Unrestricted'                
            else ''"/>
        <xsl:sequence select="$methodSpinType"/>        
    </xsl:function>    
    
</xsl:stylesheet>
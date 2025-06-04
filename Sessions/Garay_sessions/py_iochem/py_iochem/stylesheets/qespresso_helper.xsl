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
    xmlns:qex="http://www.iochem-bd.org/dictionary/qespresso/"
    xmlns:helper="http://www.w3.org/1999/XSL/Helper-Functions"
    xmlns:cml="http://www.xml-cml.org/schema"
    xmlns:cmlx="http://www.xml-cml.org/schema/cmlx"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    exclude-result-prefixes="xs qex cml cmlx"
    version="2.0">

    <!-- Calculation type related variables, not all used but processed -->
    <xsl:variable name="qex:GeometryOptimization" select="'Geometry optimization'" />
    <xsl:variable name="qex:SinglePoint" select="'Single point'" />
    <xsl:variable name="qex:MolecularDynamics" select="'Ab-Initio Molecular Dynamics'"/>
    <xsl:variable name="qex:Bands" select="'Band structure'"/>
    <xsl:variable name="qex:NonSCF" select="'Non-SCF calculation'"/>
    <xsl:variable name="qex:CarParrinello" select="'Car-parrinello Molecular Dynamics'"/>
    <xsl:variable name="qex:CarParrinelloWF" select="'Car-parrinello Molecular Dynamics (Wannier Functions)'"/>
    <xsl:variable name="qex:NudgedElasticBand" select="'Nudged Elastic Band'"/>


    <xsl:function name="qex:getShellType">
        <xsl:param name="nspin"/>
        <xsl:variable name="nspinN" select="round(number($nspin))" />
        <xsl:choose>
            <xsl:when test="$nspinN = 1"><xsl:value-of select="'Closed'"/></xsl:when>
            <xsl:when test="$nspinN = 2"><xsl:value-of select="'Open'"/></xsl:when>
            <xsl:when test="$nspinN = 4"><xsl:value-of select="'Open (Noncollinear)'"/></xsl:when>
            <xsl:otherwise><xsl:value-of select="'Closed'"/></xsl:otherwise>
        </xsl:choose>
    </xsl:function>

    <xsl:function name="qex:getBasis">
        <xsl:param name="pseudos" as="xs:string*"/>

        <!-- Basis abbreviations -->
        <xsl:variable name="qex:pseudos">
            <entry key="ULTRASOFT" value="Ultrasoft"/>
            <entry key="NORM-CONSERVING" value="NORMCONS"/>
            <entry key="PROJECTOR AUGMENTED-WAVE" value="PAW"/>
            <entry key="COULOMB" value="Coulomb Potential"/>
        </xsl:variable>
        <xsl:for-each select="$pseudos">
            <xsl:variable name="pseudo" select="."/>
            <xsl:value-of select="$qex:pseudos/*:entry[contains(upper-case($pseudo),@key)]/@value"/>
        </xsl:for-each>
    </xsl:function>

    <xsl:function name="qex:getCharge">
        <xsl:param name="charge"/>
        <xsl:choose>
            <xsl:when test="not(exists($charge)) or $charge = ''">0</xsl:when>
            <xsl:otherwise><xsl:value-of select="number($charge)"/></xsl:otherwise>
        </xsl:choose>
    </xsl:function>

    <xsl:function name="qex:getCalcType">
        <xsl:param name="moduleName"/>
        <xsl:param name="calculation"/>

        <xsl:variable name="nCalculation" select="replace(helper:trim(upper-case($calculation)),'[^A-Z-]','')"/>

        <xsl:choose>
            <xsl:when test="helper:trim(upper-case($moduleName)) = 'PWSCF'">
                <xsl:choose>
                    <xsl:when test="$nCalculation = 'VC-RELAX'">
                        <xsl:value-of select="$qex:GeometryOptimization"/>
                    </xsl:when>
                    <xsl:when test="$nCalculation = 'RELAX'">
                        <xsl:value-of select="$qex:GeometryOptimization"/>
                    </xsl:when>
                    <xsl:when test="$nCalculation = 'SCF'">
                        <xsl:value-of select="$qex:SinglePoint"/>
                    </xsl:when>
                    <xsl:when test="$nCalculation = 'BANDS'">
                        <xsl:value-of select="$qex:Bands"/>
                    </xsl:when>
                    <xsl:when test="$nCalculation = 'NSCF'">
                        <xsl:value-of select="$qex:NonSCF"/>
                    </xsl:when>
                    <xsl:when test="$nCalculation = 'MD'">
                        <xsl:value-of select="$qex:MolecularDynamics"/>
                    </xsl:when>
                    <xsl:when test="$nCalculation = 'CP'">
                        <xsl:value-of select="$qex:CarParrinello"/>
                    </xsl:when>
                    <xsl:when test="$nCalculation = 'CP-WF'">
                        <xsl:value-of select="$qex:CarParrinelloWF"/>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:value-of select="$qex:SinglePoint"/>
                    </xsl:otherwise>
                </xsl:choose>
            </xsl:when>
            <xsl:when test="helper:trim(upper-case($moduleName)) = 'PWNEB'">
                <xsl:value-of select="$qex:NudgedElasticBand"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="$qex:SinglePoint"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:function>

   <!-- Do no change the order of the following terms inside regex, it can make the functions fail -->
    <xsl:variable name="functionalRegex" select="'.*(star1s-blyp|star1s-pbe|star1s-pz|blyp|bp|coulomb|rel-pbesol|rel-pbe|pbesol|pbe|pw91|rel-pz|pz|tpss).*'"/>
    <xsl:variable name="methodRegex" select="'.*(ae|bhs|bpaw|hgh|kjpaw|mt|nrrkjus|rrkjus|rrkj|van|vbc).*'"/>

    <xsl:function name="qex:getFunctionalsFromFilenames">
        <xsl:param name="filenames" as="xs:string*"/>
        <xsl:for-each select="$filenames">
            <xsl:analyze-string select="lower-case(.)" regex="{$functionalRegex}">
                <xsl:matching-substring>
                    <xsl:value-of select="regex-group(1)"/>
                </xsl:matching-substring>
            </xsl:analyze-string>
        </xsl:for-each>
    </xsl:function>

    <xsl:function name="qex:getMethodsFromFilenames">
        <xsl:param name="filenames" as="xs:string*"/>
        <xsl:for-each select="$filenames">
            <xsl:analyze-string select="lower-case(.)" regex="{$methodRegex}">
                <xsl:matching-substring>
                    <xsl:value-of select="regex-group(1)"/>
                </xsl:matching-substring>
            </xsl:analyze-string>
        </xsl:for-each>
    </xsl:function>

    <xsl:function name="qex:getInputParameter">
        <xsl:param name="path"/>
        <xsl:param name="parameterName"/>
        <xsl:value-of select="$path//cml:module[@dictRef='cc:initialization']//cml:parameter[@dictRef='input']//cml:scalar[@dictRef='cc:parameter' and text() = $parameterName]/following-sibling::cml:scalar[@dictRef='cc:value']"/>
    </xsl:function>


    <xsl:function name="qex:getInputParameter">
        <xsl:param name="path"/>
        <xsl:param name="templateRef"/>
        <xsl:param name="parameterName"/>
        <xsl:value-of select="$path//cml:module[@dictRef='cc:initialization']//cml:parameter[@dictRef='input']/cml:list[upper-case(@cmlx:templateRef)= upper-case($templateRef)]//cml:scalar[@dictRef='cc:parameter' and text() = $parameterName]/following-sibling::cml:scalar[@dictRef='cc:value']"/>
    </xsl:function>

</xsl:stylesheet>

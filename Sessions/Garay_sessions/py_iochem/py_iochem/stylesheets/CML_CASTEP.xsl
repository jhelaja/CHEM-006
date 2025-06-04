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
<xsl:stylesheet
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xd="http://www.oxygenxml.com/ns/doc/xsl"
    xmlns:cml="http://www.xml-cml.org/schema"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:cmlx="http://www.xml-cml.org/schema/cmlx"
    xmlns:ckbk="http://my.safaribooksonline.com/book/xml/0596009747/numbers-and-math/77"
    xmlns:ca="http://www.iochem-bd.org/dictionary/castep/"
    xmlns:helper="http://www.w3.org/1999/XSL/Helper-Functions"

    xpath-default-namespace="http://www.xml-cml.org/schema" exclude-result-prefixes="xs xd cml ckbk ca helper cmlx"
    version="2.0">
    <xd:doc scope="stylesheet">
        <xd:desc>
            <xd:p><xd:b>Created on:</xd:b> May 29 2024</xd:p>
            <xd:p><xd:b>Author:</xd:b>Moisés Álvarez Moreno and Diego Garay Ruiz</xd:p>
            <xd:p><xd:b>Center:</xd:b>Institute of Chemical Research of Catalonia</xd:p>
        </xd:desc>
    </xd:doc>
    <xsl:include href="chem_helper.xsl"/>
    <xsl:include href="castep_helper.xsl"/>
    <xsl:output method="text" omit-xml-declaration="yes" indent="no" />
    <xsl:strip-space elements="*"/>

    <xsl:param name="title"/>
    <xsl:param name="author"/>
    
    <!-- Environment module -->
    <xsl:variable name="environment" select="//cml:module[@id='job'][1]/cml:module[@id='environment']/cml:parameterList"/>
    <xsl:variable name="programParameter" select="$environment/cml:parameter[@dictRef='cc:program']"/>
    <xsl:variable name="versionParameter" select="$environment/cml:parameter[@dictRef='cc:programVersion']"/>

    <!-- Geometry -->
    <xsl:variable name="initialMolecule" select="(//cml:molecule[@id='initial'])[last()]"/>
    <xsl:variable name="finalMolecule" select="(//cml:molecule[@id='final'])[last()]"/>

    <!-- Initialization -->
    <xsl:variable name="initialization" select="//cml:module[@id='job'][1]/cml:module[@dictRef='cc:initialization']"/>
    <xsl:variable name="setup" select="$initialization/cml:parameterList" />
    <xsl:variable name="calcType" select="ca:getCalcType($setup)" />
    <xsl:variable name="isOptimization" select="contains($calcType,$ca:minimize)"/>
    <xsl:variable name="pointGroup" select="$initialization/cml:module[@id='otherComponents']/cml:module[@cmlx:templateRef='symmetry']" />
    <xsl:variable name="kpoints" select="$initialization/cml:module[@id='otherComponents']/cml:module[@cmlx:templateRef='kpoints']/cml:matrix" />

    <!-- Calculation -->
    <xsl:variable name="calculation" select="//cml:module[@id='job'][1]/cml:module[@dictRef='cc:calculation']"/>

    <!-- Finalization -->
    <xsl:variable name="finalization" select="//cml:module[@id='job'][1]/cml:module[@dictRef='cc:finalization']"/>
    <xsl:variable name="graphs" select="$finalization/cml:module[@dictRef='cc:userDefinedModule']/cml:module[@id='chart']" />

    <xsl:template match="/">
 
        <xsl:call-template name="generalInfo"/>
        <xsl:call-template name="settings"/>
        <xsl:call-template name="atomicCoordinates" />
                
        <xsl:if test="exists($kpoints)">
                    <xsl:call-template name="pointGroup"/>
                    <xsl:call-template name="kpoints"/>
        </xsl:if>

        <xsl:for-each select="//cml:module[@dictRef='cc:job']">
            <xsl:call-template name="energies" />
            <xsl:call-template name="timing"/>
            <xsl:text> FETCH_JOB_END&#xa;</xsl:text>
<!--                    <xsl:call-template name="mulliken" />-->
                        
        </xsl:for-each>
                    
    </xsl:template>

    <!-- Atomic coordinates -->
    <xsl:template name="atomicCoordinates">
        <xsl:variable name="isLargeMolecule" select="count($finalMolecule/atomArray/atom) > 1000"/>
        <xsl:variable name="collapseAccordion" select="if($isLargeMolecule) then '' else 'in'"/>

        <xsl:if test="$isOptimization">
            <xsl:choose>
                <xsl:when test="exists($calculation//cml:module[@cmlx:templateRef='scf']//cml:scalar[@dictRef='cc:scfConverged']) and not(exists(//cml:molecule[@id='not.converged']))">
                    <xsl-text>optimizationState:optimized#;#</xsl-text>
                </xsl:when>
                <xsl:otherwise>
                    <xsl-text>optimizationState:not_converged#;#</xsl-text>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:if>
              
        <xsl-text> geometryCartesian: </xsl-text>
        <xsl:for-each select="$finalMolecule/cml:atomArray/cml:atom">
            <xsl:variable name="elementType" select="./@elementType"/>
            <xsl:value-of select="concat($elementType,' ',@x3,' ',@y3,' ',@z3,'&#xa;')"/>
        </xsl:for-each>
        <xsl-text> #;# </xsl-text>
    </xsl:template>

    <!-- General Info -->
    <xsl:template name="generalInfo">

        <xsl:if test="$title">
            <xsl:value-of select="$title"/>
        </xsl:if>
        <xsl:value-of select="concat('program:',$programParameter,'#;#')"/>
        <xsl:value-of select="concat('version:',$versionParameter,'#;#')"/>
        <xsl:if test="$author">
            <xsl:value-of select="concat('author:',$author,'#;#')"/>
        </xsl:if>
        <xsl:value-of select="concat('formula:',(//formula/@concise)[1],'#;#')"/>
    </xsl:template>

    <!-- Settings -->
    <xsl:template name="settings">
            <xsl:for-each select="('Geometry Optimization', 'Electronic Spectroscopy', 'Exchange Correlation')">
                <xsl:variable name="header" select="." />
                <xsl:variable name="section" select="replace(lower-case(.),' ','.')"/>
                <xsl:variable name="tablename" select="replace(lower-case(.), ' ', '')" />
                <xsl:if test="exists($setup/cml:parameter[@title=$section])">
                    <xsl:copy-of select="helper:printSectionRows($setup, $section)"/>
                </xsl:if>
            </xsl:for-each>
    </xsl:template>

    <!-- Kpints -->
    <xsl:template name="kpoints">
        <xsl:if test="exists($kpoints)">
            <xsl:variable name="kpmatrix" select="tokenize($kpoints/text(), '\s+')"/>
            <xsl-text> kpoints: </xsl-text>
            <xsl:for-each select="1 to $kpoints/@rows">
                <xsl:variable name="outerIndex" select="."/>
                <xsl:variable name="kp1" select="number($kpmatrix[(($outerIndex - 1) * $kpoints/@rows) + 1])"/>
                <xsl:variable name="kp2" select="number($kpmatrix[(($outerIndex - 1) * $kpoints/@rows) + 1])"/>
                <xsl:variable name="kp3" select="number($kpmatrix[(($outerIndex - 1) * $kpoints/@rows) + 1])"/>
                <xsl:variable name="kp4" select="number($kpmatrix[(($outerIndex - 1) * $kpoints/@rows) + 1])"/>
                <xsl:value-of select="concat($kp1,' ',$kp2,' ',$kp3,' ',$kp4,'&#xa;')"/>
            </xsl:for-each> 
        </xsl:if>
    </xsl:template>

    <!-- Point group line -->
    <xsl:template name="pointGroup">
        <xsl:if test="exists($pointGroup)">
            <xsl:value-of select="concat('pointGroup:',$pointGroup//cml:scalar[@dictRef='cc:pointgroup'],'#;#')"/>
            <xsl:value-of select="concat('spaceGroup:',$pointGroup//cml:scalar[@dictRef='ca:spacegroup'],'#;#')"/>
        </xsl:if>
    </xsl:template>

    <!-- Modules -->
    <xsl:template name="energies">
        <xsl:variable name="scf" select="($calculation/cml:module[@id='otherComponents']/cml:module[@cmlx:templateRef='step']/cml:module[@cmlx:templateRef='scf'])[last()]" />
        <xsl:if test="exists($scf)">
            <xsl:variable name="finalEnergy" select="$scf/cml:scalar[@dictRef='cc:finalEnergy']" />
            <xsl:variable name="freeEnergy" select="$scf/cml:scalar[@dictRef='cc:freeEnergy']" />
            <xsl:variable name="enthalpy" select="$finalization/cml:propertyList/cml:property[@dictRef='ca:enthalpy']/cml:scalar" />
            <xsl:variable name="converged" select="$scf/cml:scalar[@dictRef='cc:scfConverged']" />
            
            <xsl:value-of select="concat('finalEnergy:',$finalEnergy,'#;#')"/>
            <xsl:value-of select="concat('finalEnergyUnits:',helper:printUnitSymbol($finalEnergy/@units),'#;#')"/>
            <xsl:value-of select="concat('freeEnergy:',$freeEnergy,'#;#')"/>
            <xsl:value-of select="concat('freeEnergyUnits:',helper:printUnitSymbol($freeEnergy/@units),'#;#')"/>
            <xsl:value-of select="concat('enthalpy:',$enthalpy,'#;#')"/>
            <xsl:value-of select="concat('enthalpyUnits:',helper:printUnitSymbol($enthalpy/@units),'#;#')"/>
            <xsl:value-of select="concat('scfConvergence:',$converged,'#;#')"/>
            
        </xsl:if>
    </xsl:template>

<!--    <xsl:template name="mulliken">
        <xsl:variable name="mulliken" select="$finalization/cml:module[@id='otherComponents']/cml:module[@cmlx:templateRef='mulliken']" />
        <xsl:if test="exists($mulliken)">

            <xsl:variable name="serial" select="tokenize($mulliken/list/*[@dictRef='cc:serial'],'\s+')"/>
            <xsl:variable name="elementType" select="tokenize($mulliken/list/*[@dictRef='cc:elementType'],'\s+')"/>
            <xsl:variable name="spinType" select="tokenize($mulliken/list/*[@dictRef='ca:spin'],'\s+')"/>
            <xsl:variable name="orbitalS" select="tokenize($mulliken/list/*[@dictRef='ca:orbitalS'],'\s+')"/>
            <xsl:variable name="orbitalP" select="tokenize($mulliken/list/*[@dictRef='ca:orbitalP'],'\s+')"/>
            <xsl:variable name="orbitalD" select="tokenize($mulliken/list/*[@dictRef='ca:orbitalD'],'\s+')"/>
            <xsl:variable name="orbitalF" select="tokenize($mulliken/list/*[@dictRef='ca:orbitalF'],'\s+')"/>
            <xsl:variable name="total" select="tokenize($mulliken/list/*[@dictRef='x:total'],'\s+')"/>
            <xsl:variable name="charge" select="tokenize($mulliken/list/*[@dictRef='x:charge'],'\s+')"/>
            <xsl:variable name="spin" select="tokenize($mulliken/list/*[@dictRef='x:spin'],'\s+')"/>

            <div class="panel panel-default">
                <div class="panel-heading" data-toggle="collapse" data-target="div#charges-{generate-id($mulliken)}" style="cursor: pointer; cursor: hand;">
                    <h4 class="panel-title">
                        Mulliken Atomic Charges
                    </h4>
                </div>
                <div id="charges-{generate-id($mulliken)}" class="panel-collapse collapse">
                    <div class="panel-body">
                        <div class="row bottom-buffer">
                            <div class="col-lg-8 col-md-8 col-sm-12">
                                <table id="chargesT-{generate-id($mulliken)}"></table>
                                <script type="text/javascript">
                                    $(document).ready(function() {
                                        $('table#chargesT-<xsl:value-of select="generate-id($mulliken)"/>').dataTable({
                                            "aaData": [
                                            <xsl:for-each select="1 to count($elementType)">
                                                <xsl:variable name="innerIndex" select="."/>
                                                [ <xsl:value-of select="$serial[$innerIndex * 2 - 1]"/>,
                                                  "<xsl:value-of select="$elementType[$innerIndex]"/>",
                                                  "<xsl:value-of select="$spinType[$innerIndex * 2 - 1]"/>",
                                                  "<xsl:value-of select="$orbitalS[$innerIndex * 2 - 1]"/>",
                                                  "<xsl:value-of select="$orbitalP[$innerIndex * 2 - 1]"/>",
                                                  "<xsl:value-of select="$orbitalD[$innerIndex * 2 - 1]"/>",
                                                  "<xsl:value-of select="$orbitalF[$innerIndex * 2 - 1]"/>",
                                                  "<xsl:value-of select="$total[$innerIndex * 2 - 1]"/>",
                                                  "<xsl:value-of select="$charge[$innerIndex]"/>",
                                                  "<xsl:value-of select="$spin[$innerIndex]"/>"
                                                    ],
                                                [ "<xsl:value-of select="$serial[$innerIndex * 2]"/>",
                                                  '',
                                                  "<xsl:value-of select="$spinType[$innerIndex * 2]"/>",
                                                  "<xsl:value-of select="$orbitalS[$innerIndex * 2]"/>",
                                                  "<xsl:value-of select="$orbitalP[$innerIndex * 2]"/>",
                                                  "<xsl:value-of select="$orbitalD[$innerIndex * 2]"/>",
                                                  "<xsl:value-of select="$orbitalF[$innerIndex * 2]"/>",
                                                  "<xsl:value-of select="$total[$innerIndex * 2 - 1]"/>",
                                                  '',
                                                  '']
                                                <xsl:if test="($innerIndex &lt; count($elementType))"><xsl:text>,</xsl:text></xsl:if>

                                            </xsl:for-each>
                                            ],
                                            "aoColumns": [
                                                { "sTitle": "Atom" },
                                                { "sTitle": ""},
                                                { "sTitle": ""},
                                                { "sTitle": "S" },
                                                { "sTitle": "P" },
                                                { "sTitle": "D" },
                                                { "sTitle": "F"  },
                                                { "sTitle": "Total"  },
                                                { "sTitle": "Charge"},
                                                { "sTitle": "Spin"  },
                                            ],
                                            "aoColumnDefs" : [
                                                { "sClass": "text-right", "aTargets": [3,4,5,6,7,8,9] }
                                            ],
                                            "bFilter": false,
                                            "bPaginate": true,
                                            "bSort": false,
                                            "bInfo": true
                                            } );
                                    } );
                                </script>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </xsl:if>

    </xsl:template>-->


    <!-- Timing -->
    <xsl:template name="timing">
        <xsl:variable name="inittime" select="$finalization/cml:propertyList/cml:property[@dictRef='ca:inittime']" />
        <xsl:variable name="calctime" select="$finalization/cml:propertyList/cml:property[@dictRef='ca:calctime']" />
        <xsl:variable name="endtime" select="$finalization/cml:propertyList/cml:property[@dictRef='ca:endtime']" />
        <xsl:variable name="wallTime" select="$finalization/cml:propertyList/cml:property[@dictRef='cc:wallTtime']" />
      
        <xsl:if test="exists($wallTime)">
            <xsl:value-of select="concat('totalTime:',$walltime,'#;#')"/>
            <xsl:value-of select="concat('totalTimeUnits:',helper:printUnitSymbol($wallTime/cml:scalar/@units),'#;#')"/>
        </xsl:if>
    </xsl:template>

    <xsl:function name="helper:printSectionRows">
        <xsl:param name="params" />
        <xsl:param name="title" />

        <xsl:for-each select="$params/cml:parameter[@title=$title]">
            <xsl:variable name="param" select="."/>
            <xsl:variable name="paramLabel" select="replace($param/cml:scalar[@dictRef='x:label'],' ','_')"/>
            <xsl:variable name="paramValue" select="replace($param/cml:scalar[@dictRef='x:value'],' ','_')"/>
            <xsl:value-of select="concat($paramLabel,':',$paramValue,'#;#')"/>
            <xsl:if test="exists($param/cml:scalar[@dictRef='x:value']/@units)">
                <xsl:value-of select="concat($paramLabel,'Units:',helper:printUnitSymbol($param/cml:scalar[@dictRef='x:value']/@units),'#;#')"/>
            </xsl:if>
        </xsl:for-each>
    </xsl:function>

</xsl:stylesheet>

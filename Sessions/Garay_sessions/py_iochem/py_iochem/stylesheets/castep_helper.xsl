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
    xmlns:ca="http://www.iochem-bd.org/dictionary/castep/"
    xmlns:cml="http://www.xml-cml.org/schema"
    xmlns:cmlx="http://www.xml-cml.org/schema/cmlx"
    xmlns:helper="http://www.w3.org/1999/XSL/Helper-Functions"
    exclude-result-prefixes="xs"
    version="2.0">

    <xsl:variable name="ca:minimize">Geometry optimization</xsl:variable>
    <xsl:variable name="ca:spectroscopy">Electronic spectroscopy</xsl:variable>

    <xsl:function name="ca:getCalcType">
        <xsl:param name="setup"  />

        <xsl:variable name="type" select="ca:getParameter($setup, 'type of calculation')" />
        <xsl:choose>
            <xsl:when test="matches($type, 'geometry optimization')">
                <xsl:value-of select="$ca:minimize"/>
            </xsl:when>
            <xsl:when test="matches($type, 'Electronic Spectroscopy')">
                <xsl:value-of select="$ca:spectroscopy"/>
            </xsl:when>
            <xsl:otherwise><xsl:text>N/A</xsl:text></xsl:otherwise>
        </xsl:choose>
    </xsl:function>

    <xsl:function name="ca:getParameter">
        <xsl:param name="setup" />
        <xsl:param name="name" />
        <xsl:copy-of select="$setup/cml:parameter/cml:scalar[@dictRef='x:label'][text()=$name]/following-sibling::cml:scalar[@dictRef='x:value']"/>
    </xsl:function>

</xsl:stylesheet>

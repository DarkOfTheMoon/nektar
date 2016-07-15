<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Tetrahedral elements, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_channel_robin_m3.xml</parameters>
    <files>
        <file description="Session File">Tet_channel_robin_m3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000788723</value>
            <value variable="v" tolerance="1e-12">8.39367e-05</value>
            <value variable="w" tolerance="1e-12">0.000450986</value>
	    <value variable="p" tolerance="1e-12">0.0156638</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00344169</value>
            <value variable="v" tolerance="1e-12">0.000306357</value>
            <value variable="w" tolerance="1e-12">0.0018145</value>
	    <value variable="p" tolerance="1e-12">0.0866085</value>
        </metric>
    </metrics>
</test>


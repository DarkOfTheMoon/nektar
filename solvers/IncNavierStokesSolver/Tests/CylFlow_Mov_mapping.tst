<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving cylinder, 2D, using Mappings for forced oscillations</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_Mov_mapping.xml</parameters>
    <files>
        <file description="Session File">CylFlow_Mov_mapping.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">21.4188</value>
            <value variable="v" tolerance="1e-9">0.956901</value>
	    <value variable="p" tolerance="1e-9">2.43895</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">1.88686</value>
            <value variable="v" tolerance="1e-9">0.911823</value>
	    <value variable="p" tolerance="1e-9">1.75362</value>
        </metric>
    </metrics>
</test>



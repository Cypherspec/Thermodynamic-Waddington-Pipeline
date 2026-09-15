# CAD assembly BOM — Revision B

Quantities are installed CAD occurrences, not a procurement or fabrication release. All materials and vendor selections are pending. Envelope parts may represent purchased assemblies with unmodeled internal components. `geometry.kcl` and `fluidRouting.kcl` are shared construction libraries, not BOM items.

| Part file | Qty | CAD role |
|---|---:|---|
| basePlate.kcl | 1 | Foundation plate |
| upright.kcl | 4 | Hollow corner posts |
| roof.kcl | 1 | Removable roof |
| sidePanel.kcl | 2 | Windowed side panels |
| sideWindow.kcl | 2 | Side glazing |
| frontDoor.kcl | 2 | Access door frames |
| doorWindow.kcl | 2 | Door glazing |
| rearPanel.kcl | 1 | Utility and fan openings |
| isolationFoot.kcl | 4 | Isolation foot envelopes |
| breadboard.kcl | 1 | Optical mounting plate |
| microscopeColumn.kcl | 1 | Microscope support |
| stageBridge.kcl | 2 | Lower and upper support bridges |
| stageRail.kcl | 2 | Guide rail envelopes |
| stagePlaten.kcl | 1 | Apertured stage |
| heaterPlate.kcl | 1 | Thermal spreader |
| dockingPlate.kcl | 1 | Dual cartridge mounting plate |
| cartridgeCarrier.kcl | 2 | Stepped slide carriers |
| coverglass.kcl | 2 | Optical cartridge bottoms |
| trapChip.kcl | 1 | Live imaging chip candidate |
| dropletChip.kcl | 1 | Encapsulation chip candidate |
| cartridgeClamp.kcl | 2 | Chip retention frames |
| clampScrew.kcl | 8 | Retention screw envelopes |
| chipFerrule.kcl | 7 | Bored tubing adapters |
| chipTube.kcl | 7 | Short tubing tails |
| incubatorHood.kcl | 1 | Environmental chamber frame |
| incubatorLid.kcl | 1 | Removable chamber lid |
| scientificCamera.kcl | 1 | Camera installation envelope |
| cameraCoupler.kcl | 1 | Bored optical spacer |
| filterCube.kcl | 1 | Bored filter enclosure envelope |
| objective.kcl | 1 | Objective barrel envelope |
| opticalElement.kcl | 2 | Flat optical window representatives |
| illuminator.kcl | 1 | Condenser housing envelope |
| illuminationArm.kcl | 1 | Condenser support |
| stageMotor.kcl | 1 | Drive installation envelope |
| lightEngine.kcl | 1 | Multicolour source envelope |
| lightGuide.kcl | 1 | Excitation coupling sleeve |
| opticalSupport.kcl | 4 | Breadboard spacers |
| spillTray.kcl | 1 | Open secondary containment tray |
| reagentRack.kcl | 1 | Six-reservoir rack |
| rackPost.kcl | 4 | Rack spacers |
| reagentBottle.kcl | 6 | Hollow reservoir envelopes |
| bottleCap.kcl | 6 | Two-port reservoir cap envelopes |
| wasteBottle.kcl | 2 | Sample and waste vessel envelopes |
| wasteLid.kcl | 2 | Collection vessel lids |
| equipmentShelf.kcl | 1 | Aft equipment shelf |
| shelfPost.kcl | 4 | Equipment shelf spacers |
| pressureController.kcl | 1 | Eight-channel controller envelope |
| valveRail.kcl | 1 | Valve/sensor shelf |
| solenoidValve.kcl | 8 | Valve envelopes |
| flowSensor.kcl | 8 | Sensor envelopes |
| computerCase.kcl | 1 | Acquisition PC enclosure |
| computerLid.kcl | 1 | PC lid |
| controllerBoard.kcl | 1 | PCB envelope |
| heatSink.kcl | 1 | Finned heat sink |
| powerSupply.kcl | 1 | Supply installation envelope |
| gasMixer.kcl | 1 | Environmental controller envelope |
| doorHandle.kcl | 2 | Curved door pulls |
| doorHinge.kcl | 4 | Hinge installation envelopes |
| displayBezel.kcl | 1 | Operator display housing |
| displayScreen.kcl | 1 | Screen insert |
| stopButton.kcl | 1 | Stop actuator envelope |
| returnHose.kcl | 1 | Hollow semicircular service-loop representative |
| cameraRiser.kcl | 1 | Camera mounting pad |
| lightPedestal.kcl | 1 | Source support bracket |
| valvePost.kcl | 2 | Valve shelf pillars |
| splashPartition.kcl | 1 | Wet/dry divider |
| exhaustFan.kcl | 1 | Fan housing envelope |

## Revision B additions

| Part file | Qty | CAD role |
|---|---:|---|
| chamberGasket.kcl | 2 | Nominal upper/lower chamber frame seals |
| chamberWindow.kcl | 1 | Separate chamber optical window |
| windowSeal.kcl | 1 | Nominal annular window seal |
| windowRetainer.kcl | 1 | Window retaining ring; attachment engineering pending |
| tubeGrommet.kcl | 7 | Chamber lid tube feedthrough inserts |
| boardStandoff.kcl | 4 | PCB mounting spacers |
| cartridgeNut.kcl | 8 | Backing nut envelopes; cosmetic thread intent |
| liveSupplyLine.kcl | 1 | Selector common to imaging cartridge |
| liveWasteLine.kcl | 1 | Imaging cartridge to waste vessel |
| sampleCollectionLine.kcl | 1 | Droplet outlet to sample vessel |
| dropletSupplyLine.kcl | 1 | Cell reservoir to droplet inlet |
| frontOilLine.kcl | 1 | Front oil reservoir to first oil inlet |
| rearOilLine.kcl | 1 | Rear oil reservoir to second oil inlet |
| beadSupplyLine.kcl | 1 | Bead reservoir to droplet chip |
| mediaSupplyLine.kcl | 1 | Media reservoir to selector input |
| perturbationLine.kcl | 1 | Perturbation reservoir to selector input |
| reagentSelector.kcl | 1 | Six-to-one selector installation envelope, not modeled internal switching |
| selectorSupport.kcl | 1 | Selector pedestal |
| reservoirGrommet.kcl | 6 | Reservoir liquid-tube port inserts |
| collectionGrommet.kcl | 2 | Sample/waste liquid-tube port inserts |
| selectorPlug.kcl | 4 | Unassigned selector input plugs |

Existing trap chip gains three underside registration recesses; incubator lid gains seven individual feedthroughs and a window counterbore; clamp screws are lengthened to engage backing nuts; illumination arm gains a liquid-feed clearance hole. Original part definitions and quantities above are otherwise retained.

## Items not yet released or fully modeled
Pressure supply and relief plumbing; gas supply/exhaust; complete liquid fittings and inline sensor connections; wiring/harnesses; exact commercial optical elements and mounts; full XY/focus motion hardware; structural and electrical mounting fasteners; safety interlock circuits; library-preparation and sequencing equipment. Nominal seals are now modeled but material, compression, retention and leakage are unvalidated. See NOTES.md and VALIDATION.md.

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#include "DetectorConstruction.hh"
#include "DR_PMTSD.hh"
#include "SurfaceProperty.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4UserLimits.hh"

#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include "CreateTree.hh"
#include <algorithm>
#include <string>
#include <sstream>
#include "G4MultiUnion.hh"

using namespace CLHEP;

DetectorConstruction::DetectorConstruction(const string &configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------

  ConfigFile config(configFileName);

  //------------- Not used --------------
  config.readInto(bar_length, "bar_length");
  config.readInto(core_radius_x, "core_radius_x");
  config.readInto(core_radius_y, "core_radius_y");
  config.readInto(core_material, "core_material");
  config.readInto(core_rIndex, "core_rIndex");
  config.readInto(core_absLength, "core_absLength");
  config.readInto(gap_size_x, "gap_size_x");
  config.readInto(gap_size_y, "gap_size_y");
  config.readInto(det_size_x, "det_size_x");
  config.readInto(det_size_y, "det_size_y");
  config.readInto(depth, "depth");
  config.readInto(cryst_dist, "cryst_dist");
  config.readInto(trackerX0, "trackerX0");
  config.readInto(services_thick, "services_thick");
  config.readInto(fiber_type, "fiber_type");
  config.readInto(scinti_material, "scinti_material");
  config.readInto(Cherenc_material, "Cherenc_material");
  config.readInto(Cherenp_material, "Cherenp_material");
  config.readInto(rear_filter, "rear_filter");
  config.readInto(front_filter, "front_filter");
  config.readInto(hole_diameter, "hole_diameter");
  config.readInto(fiber_diameter, "fiber_diameter");
  config.readInto(wrapping_thick, "wrapping_thick");
  config.readInto(hcal_width, "hcal_width");
  config.readInto(hcalTile_width, "hcalTile_width");
  config.readInto(hcalAbs_1_thick, "hcalAbs_1_thick");
  config.readInto(hcalAbs_2_thick, "hcalAbs_2_thick");
  config.readInto(solenoid_thick, "solenoid_thick");
  config.readInto(hcalTile_thick, "hcalTile_thick");

  //------------- used in Single bar --------------
  config.readInto(checkOverlaps, "checkOverlaps");
  config.readInto(world_material, "world_material");
  config.readInto(gap_l, "gap_l");
  config.readInto(gap_material, "gap_material");
  config.readInto(det_l, "det_l");
  config.readInto(det_material, "det_material");
  config.readInto(ecal_det_size, "ecal_det_size");

  config.readInto(ecal_incline, "ecal_incline");
  config.readInto(ecal_xy_gap, "ecal_xy_gap");
  config.readInto(ecal_material, "ecal_material");
  config.readInto(wrap_material, "wrap_material");
  config.readInto(wrap_ref, "wrap_ref");

  config.readInto(ecal_front_length, "ecal_front_length");
  config.readInto(ecal_rear_length, "ecal_rear_length");
  config.readInto(ecal_front_face, "ecal_front_face");
  config.readInto(ecal_rear_face, "ecal_rear_face");
  config.readInto(ecal_timing_distance, "ecal_timing_distance");

  B_field_intensity = config.read<double>("B_field_intensity") * tesla;

  expHall_x = 101. * cm;
  expHall_y = 101. * cm;
  expHall_z = 301. * cm;

  B_field_IsInitialized = false;

  initializeMaterials();
  initializeSurface();

  CreateTree::Instance()->inputTrackerX0 = trackerX0;
  CreateTree::Instance()->inputServiceAlmm = services_thick;
  CreateTree::Instance()->inputTimingThick = core_radius_x * 2;
  CreateTree::Instance()->inputE1Thick = ecal_front_length;
  CreateTree::Instance()->inputE2Thick = ecal_rear_length;
  CreateTree::Instance()->inputE1Width = ecal_front_face;
  CreateTree::Instance()->inputTimingECAL_dist = ecal_timing_distance;
}

//---- ---- ---- ---- ---- ---- ---- ---- ----  ---- ---- ---- ---- ---- ----

DetectorConstruction::~DetectorConstruction()
{
  delete stepLimit;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl;

  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------

  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid *worldS = new G4Box("worldS", 0.5*expHall_x, 0.5*expHall_y, 0.5*expHall_z);
  G4LogicalVolume *worldLV = new G4LogicalVolume(worldS, WoMaterial, "worldLV", 0, 0, 0);
  G4VPhysicalVolume *worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "worldPV", 0, false, 0, checkOverlaps);

  //********************************************
  //  Crystal bar
  //********************************************

  double window_l = 0.2 * mm;
  double Case_l = 0.01 * mm;
  double wrapper_gap = 0.1 * mm;

  double pointingAngle = ecal_incline;
  double alveola_thickness = ecal_xy_gap * mm;
  double Front_Back_distance = 10 * mm;

  // crystal
  G4Box *ecalCrystalS_f = new G4Box("ecalCrystalS_f", 0.5 * ecal_front_face, 0.5 * ecal_front_face, 0.5 * ecal_front_length);
  G4Box *ecalCrystalS_r = new G4Box("ecalCrystalS_r", 0.5 * ecal_rear_face, 0.5 * ecal_rear_face, 0.5 * ecal_rear_length);
  G4LogicalVolume *ecalCrystalL_f = new G4LogicalVolume(ecalCrystalS_f, EcalMaterial, "ecalCrystalL_f", 0, 0, 0);
  G4LogicalVolume *ecalCrystalL_r = new G4LogicalVolume(ecalCrystalS_r, EcalMaterial, "ecalCrystalL_r", 0, 0, 0);

  //Detector
  G4Box *ecalGapS = new G4Box("ecalGapS", ecal_det_size * 0.5, ecal_det_size * 0.5, 0.5 * gap_l);
  G4Box *ecalDetWindowS = new G4Box("ecalDetWindowS", ecal_det_size * 0.5, ecal_det_size * 0.5, 0.5 * window_l);
  G4Box *ecalDetS = new G4Box("ecalDetS", ecal_det_size * 0.5, ecal_det_size * 0.5, 0.5 * det_l);
  G4Box *ecalDetCaseOuterS = new G4Box("ecalDetCaseOuterS", ecal_det_size * 0.5 + Case_l, ecal_det_size * 0.5 + Case_l, 0.5 * (det_l + window_l) + Case_l);
  G4Box *ecalDetCaseInnerS = new G4Box("ecalDetCaseInnerS", ecal_det_size * 0.5, ecal_det_size * 0.5, 0.5 * (det_l + window_l) + Case_l);
  G4VSolid *ecalDetCaseS = new G4SubtractionSolid("ecalDetCaseS", ecalDetCaseOuterS, ecalDetCaseInnerS, NULL, G4ThreeVector(0, 0, Case_l));
  G4LogicalVolume *ecalGapL = new G4LogicalVolume(ecalGapS, PMTGapMaterial, "ecalGapL", 0, 0, 0);
  G4LogicalVolume *ecalDetWindowL = new G4LogicalVolume(ecalDetWindowS, WindowMaterial, "ecalDetWindowL", 0, 0, 0);
  G4LogicalVolume *ecalDetL = new G4LogicalVolume(ecalDetS, DeMaterial, "ecalDetL", 0, 0, 0);
  G4LogicalVolume *ecalDetCaseL = new G4LogicalVolume(ecalDetCaseS, WindowMaterial, "ecalDetCaseL", 0, 0, 0);

  //wrapper
  G4Box *ecalWrapper_outerS_f = new G4Box("ecalWrapper_outerS_f", 0.5 * ecal_front_face + wrapper_gap + 0.1, 0.5 * ecal_front_face + wrapper_gap + 0.1, 0.5 * ecal_front_length + wrapper_gap + 0.1);
  G4Box *ecalWrapper_innerS_f = new G4Box("ecalWrapper_innerS_f", 0.5 * ecal_front_face + wrapper_gap, 0.5 * ecal_front_face + wrapper_gap, ecal_front_length + wrapper_gap);
  G4Box *ecalWrapper_outerS_r = new G4Box("ecalWrapper_outerS_r", 0.5 * ecal_rear_face + wrapper_gap + 0.1, 0.5 * ecal_rear_face + wrapper_gap + 0.1, 0.5 * ecal_rear_length + wrapper_gap + 0.1);
  G4Box *ecalWrapper_innerS_r = new G4Box("ecalWrapper_innerS_r", 0.5 * ecal_rear_face + wrapper_gap, 0.5 * ecal_rear_face + wrapper_gap, ecal_rear_length + wrapper_gap);
  G4VSolid *ecalWrapperS_f_shell = new G4SubtractionSolid("ecalWrapperS_f_shell", ecalWrapper_outerS_f, ecalWrapper_innerS_f);
  G4VSolid *ecalWrapperS_r_shell = new G4SubtractionSolid("ecalWrapperS_r_shell", ecalWrapper_outerS_r, ecalWrapper_innerS_r);
  const G4ThreeVector arbimov  = G4ThreeVector(0,0,0);//G4ThreeVector(0, -1*(ecal_front_face - ecal_det_size) * 0.5 * cos(pointingAngle*TMath::Pi()/180.0), (ecal_front_face - ecal_det_size) * 0.5 * sin(pointingAngle*TMath::Pi()/180.0) );
  const G4ThreeVector arbimov_up  = G4ThreeVector(0,0,0);//G4ThreeVector(0, -1*(ecal_front_face - ecal_det_size) * 0.5, 0);

  G4VSolid *ecalWrapperS_f_hallow1 = new G4SubtractionSolid("ecalWrapperS_f_hallow1", ecalWrapperS_f_shell, ecalDetCaseS, NULL, G4ThreeVector(0, 0, (ecal_front_length * 0.5 + gap_l + window_l + det_l * 0.5 + Case_l))+arbimov_up);
  G4VSolid *ecalWrapperS_f_hallow2 = new G4SubtractionSolid("ecalWrapperS_f_hallow2", ecalWrapperS_f_hallow1, ecalDetCaseS, NULL, G4ThreeVector(0, 0, -1 * (ecal_front_length * 0.5 + gap_l + window_l + det_l * 0.5 + Case_l))+arbimov_up);
  G4VSolid *ecalWrapperS_f_hallow3 = new G4SubtractionSolid("ecalWrapperS_f_hallow3", ecalWrapperS_f_hallow2, ecalGapS, NULL, G4ThreeVector(0, 0, -1 * (ecal_front_length * 0.5 + gap_l * 0.5))+arbimov_up);
  G4VSolid *ecalWrapperS_f_hallow4 = new G4SubtractionSolid("ecalWrapperS_f_hallow4", ecalWrapperS_f_hallow3, ecalGapS, NULL, G4ThreeVector(0, 0, (ecal_front_length * 0.5 + gap_l * 0.5))+arbimov_up);
  G4LogicalVolume *ecalWrapperL_f = new G4LogicalVolume(ecalWrapperS_f_hallow4, WrapMaterial, "ecalWrapperL_f", 0, 0, 0);
  G4LogicalVolume *ecalWrapperL_r = new G4LogicalVolume(ecalWrapperS_r_shell, WrapMaterial, "ecalWrapperL_r", 0, 0, 0);

  // ECAL physical placement
  const int NECAL_CRYST = 1; //6400; //2500;
  G4VPhysicalVolume *ecalCrystalP_f[NECAL_CRYST];
  //G4VPhysicalVolume *ecalCrystalP_r[NECAL_CRYST];
  G4VPhysicalVolume *ecalWrapperP_f[NECAL_CRYST];
  //G4VPhysicalVolume *ecalWrapperP_r[NECAL_CRYST];
  G4VPhysicalVolume *ecalGapP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalGapP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetWindowP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetWindowP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetCaseP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetCaseP_fr[NECAL_CRYST];

  G4VPhysicalVolume *ecalDetP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetP_fr[NECAL_CRYST];

  G4LogicalSkinSurface *ecalWrapperSurface_f = new G4LogicalSkinSurface("ecalWrapperSurface_f", ecalWrapperL_f, fFiberWrapSurface);
  G4LogicalSkinSurface *DetSurface = new G4LogicalSkinSurface("DetSurface", ecalDetL, fPMTSurface);
  G4LogicalSkinSurface *DetCase_Surface = new G4LogicalSkinSurface("DetCase_Surface", ecalDetCaseL, fPMTCaseSurface);
  G4LogicalSkinSurface *crystalSurface = new G4LogicalSkinSurface("crystalSurface", ecalCrystalL_f, fIdealPolishedOpSurface);
  //G4LogicalBorderSurface* FilterSurface_ff[NECAL_CRYST]; 

  char name[60];
  G4double x_pos[NECAL_CRYST];
  G4double y_pos[NECAL_CRYST];
  int nArrayECAL = (int)sqrt(NECAL_CRYST);
  int iCrystal;
  for (int iX = 0; iX < nArrayECAL; iX++)
  {
    for (int iY = 0; iY < nArrayECAL; iY++)
    {

      G4RotationMatrix *piRotEcal = new G4RotationMatrix;
      piRotEcal->rotateX(pointingAngle * deg);
      G4RotationMatrix *piRotEcal_op = new G4RotationMatrix;
      piRotEcal_op->rotateX((pointingAngle + 180) * deg);

      iCrystal = nArrayECAL * iX + iY;
      x_pos[iCrystal] = (iX - nArrayECAL / 2) * (ecal_front_face + alveola_thickness);
      y_pos[iCrystal] = (iY - nArrayECAL / 2) * (ecal_front_face / cos(pointingAngle * deg) + alveola_thickness);
      if (nArrayECAL / 2 != 0)
      {
        x_pos[iCrystal] -= alveola_thickness;
        y_pos[iCrystal] -= alveola_thickness;
      }
      cout << " x_pos [" << iCrystal << "] = " << x_pos[iCrystal] << " :: y_pos[" << iCrystal << "] = " << y_pos[iCrystal] << " :: angle = [" << pointingAngle * iX << ", " << pointingAngle * iY << "] " << endl;

      //---------------------- crystal and wrapper --------------------------
      sprintf(name, "ecalCrystalP_f_%d", iCrystal);
      ecalCrystalP_f[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length * 0.5), ecalCrystalL_f, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalWrapperP_f_%d", iCrystal);
      ecalWrapperP_f[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length * 0.5), ecalWrapperL_f, name, worldLV, false, 0, checkOverlaps);
/*
      //---------------------- detectors at both ends --------------------------
      sprintf(name, "ecalDetP_ff_%d", iCrystal);
      ecalDetP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l + window_l + det_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l + window_l + det_l * 0.5) * cos(pointingAngle * deg))+arbimov, ecalDetL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetP_fr_%d", iCrystal);
      ecalDetP_fr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_front_length + gap_l + window_l + det_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l + window_l + det_l * 0.5) * cos(pointingAngle * deg))+arbimov, ecalDetL, name, worldLV, false, 0, checkOverlaps);

      //---------------------- gaps, SiPM at both ends --------------------------
      sprintf(name, "ecalGapP_ff_%d", iCrystal);
      ecalGapP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l * 0.5) * cos(pointingAngle * deg))+arbimov, ecalGapL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalGapP_fr_%d", iCrystal);
      ecalGapP_fr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_front_length + gap_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l * 0.5) * cos(pointingAngle * deg))+arbimov, ecalGapL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetWindowP_ff_%d", iCrystal);
      ecalDetWindowP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l + window_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l + window_l * 0.5) * cos(pointingAngle * deg))+arbimov, ecalDetWindowL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetWindowP_fr_%d", iCrystal);
      ecalDetWindowP_fr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_front_length + gap_l + window_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l + window_l * 0.5) * cos(pointingAngle * deg))+arbimov, ecalDetWindowL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetCaseP_ff_%d", iCrystal);
      ecalDetCaseP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l + window_l + det_l * 0.5 + Case_l) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l + window_l + det_l * 0.5 + Case_l) * cos(pointingAngle * deg))+arbimov, ecalDetCaseL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetCaseP_fr_%d", iCrystal);
      ecalDetCaseP_fr[iCrystal] = new G4PVPlacement(piRotEcal_op, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_front_length + gap_l + window_l + det_l * 0.5 + Case_l) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l + window_l + det_l * 0.5 + Case_l) * cos(pointingAngle * deg))+arbimov, ecalDetCaseL, name, worldLV, false, 0, checkOverlaps);
*/
    }
  }

  //-----------------------------------------------------
  //------------- sensitive detector & surface --------------
  //-----------------------------------------------------

  auto sdman = G4SDManager::GetSDMpointer();
  auto fDRDetSD = new DR_PMTSD("/DR_Det");
  sdman->AddNewDetector(fDRDetSD);
  ecalDetL->SetSensitiveDetector(fDRDetSD);

  //-----------------------------------------------------
  //------------- Visualization attributes --------------
  //-----------------------------------------------------

  G4Colour white(1.00, 1.00, 1.00); // white
  G4Colour gray(0.50, 0.50, 0.50);  // gray
  G4Colour black(0.00, 0.00, 0.00); // black
  G4Colour red(1.00, 0.00, 0.00);   // red
  G4Colour green(0.00, 1.00, 0.00); // green
  G4Colour blue(0.00, 0.00, 1.00);  // blue
  G4Colour lightblue(0.678, 0.847, 0.902);
  G4Colour cyan(0.00, 1.00, 1.00);    // cyan
  G4Colour air(0.90, 0.94, 1.00);     // cyan
  G4Colour magenta(1.00, 0.00, 1.00); // magenta
  G4Colour yellow(1.00, 1.00, 0.00);  // yellow
  G4Colour brass(0.80, 0.60, 0.40);   // brass
  G4Colour brown(0.70, 0.40, 0.10);   // brown

  G4VisAttributes *VisAttWorld = new G4VisAttributes(white);
  VisAttWorld->SetVisibility(true);
  VisAttWorld->SetForceWireframe(true);
  worldLV->SetVisAttributes(VisAttWorld);

  G4VisAttributes *VisAttCore = new G4VisAttributes(green);
  VisAttCore->SetVisibility(false);
  VisAttCore->SetForceWireframe(true);

  G4VisAttributes *VisCrystalCore = new G4VisAttributes(brass);
  VisCrystalCore->SetVisibility(false);
  VisCrystalCore->SetForceWireframe(true);

  G4VisAttributes *VisAttGap = new G4VisAttributes(blue);
  VisAttGap->SetVisibility(false);
  VisAttGap->SetForceWireframe(true);

  G4VisAttributes *VisAttecalDet = new G4VisAttributes(gray);
  VisAttecalDet->SetVisibility(true);
  VisAttecalDet->SetForceWireframe(true);
  ecalDetL->SetVisAttributes(VisAttecalDet);
  ecalGapL->SetVisAttributes(VisAttecalDet);

  G4VisAttributes *VisAttecalDetW = new G4VisAttributes(red);
  VisAttecalDetW->SetVisibility(true);
  VisAttecalDetW->SetForceWireframe(true);
  ecalDetWindowL->SetVisAttributes(VisAttecalDetW);

  G4VisAttributes *Visabs = new G4VisAttributes(brass);
  Visabs->SetVisibility(true);
  Visabs->SetForceWireframe(true);

  G4VisAttributes *Viswrap = new G4VisAttributes(cyan);
  Viswrap->SetVisibility(true);
  Viswrap->SetForceWireframe(true);
  ecalWrapperL_f->SetVisAttributes(Viswrap);
  ecalWrapperL_r->SetVisAttributes(Viswrap);

  G4VisAttributes *VisFiber = new G4VisAttributes(lightblue);
  VisFiber->SetVisibility(true);
  VisFiber->SetForceWireframe(false);
  ecalCrystalL_f->SetVisAttributes(VisFiber);
  ecalCrystalL_r->SetVisAttributes(VisFiber);

  if (B_field_intensity > 0.1 * tesla)
    ConstructField();

  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl;
  return worldPV;
}

void DetectorConstruction::initializeMaterials()
{
  //-----------------
  // define materials

  WoMaterial = NULL;
  if (world_material == 0)
    WoMaterial = MyMaterials::Vacuum();
  else if (world_material == 1)
    WoMaterial = MyMaterials::Air();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl;
    exit(-1);
  }
  G4cout << "Wo. material: " << WoMaterial << G4endl;

  CoMaterial = NULL;
  if (core_material == 1)
    CoMaterial = MyMaterials::Quartz();
  else if (core_material == 2)
    CoMaterial = MyMaterials::SiO2();
  else if (core_material == 3)
    CoMaterial = MyMaterials::SiO2_Ce();
  else if (core_material == 4)
    CoMaterial = MyMaterials::LuAG_Ce();
  else if (core_material == 5)
    CoMaterial = MyMaterials::YAG_Ce();
  else if (core_material == 6)
    CoMaterial = MyMaterials::LSO();
  else if (core_material == 7)
    CoMaterial = MyMaterials::LYSO();
  else if (core_material == 8)
    CoMaterial = MyMaterials::LuAG_undoped();
  else if (core_material == 9)
    CoMaterial = MyMaterials::GAGG_Ce();
  else if (core_material == 11)
    CoMaterial = MyMaterials::LuAG_Pr();
  else if (core_material == 12)
    CoMaterial = MyMaterials::PbF2();
  else if (core_material == 13)
    CoMaterial = MyMaterials::PlasticBC408();
  else if (core_material == 14)
    CoMaterial = MyMaterials::PlasticBC418();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << core_material << G4endl;
    exit(-1);
  }
  G4cout << "Co. material: " << CoMaterial << G4endl;

  EcalMaterial = NULL;
  if (ecal_material == 1)
    EcalMaterial = MyMaterials::Quartz();
  else if (ecal_material == 2)
    EcalMaterial = MyMaterials::SiO2();
  else if (ecal_material == 3)
    EcalMaterial = MyMaterials::SiO2_Ce();
  else if (ecal_material == 4)
    EcalMaterial = MyMaterials::LuAG_Ce();
  else if (ecal_material == 5)
    EcalMaterial = MyMaterials::YAG_Ce();
  else if (ecal_material == 6)
    EcalMaterial = MyMaterials::LSO();
  else if (ecal_material == 7)
    EcalMaterial = MyMaterials::LYSO();
  else if (ecal_material == 8)
    EcalMaterial = MyMaterials::LuAG_undoped();
  else if (ecal_material == 9)
    EcalMaterial = MyMaterials::GAGG_Ce();
  else if (ecal_material == 10)
    EcalMaterial = MyMaterials::LuAG_Pr();
  else if (ecal_material == 11)
    EcalMaterial = MyMaterials::PbF2();
  else if (ecal_material == 12)
    EcalMaterial = MyMaterials::PlasticBC408();
  else if (ecal_material == 13)
    EcalMaterial = MyMaterials::PlasticBC418();
  else if (ecal_material == 14)
    EcalMaterial = MyMaterials::PWO();
  else if (ecal_material == 15)
    EcalMaterial = MyMaterials::Acrylic();
  else if (ecal_material == 16)
    EcalMaterial = MyMaterials::copper();
  else if (ecal_material == 17)
    EcalMaterial = MyMaterials::EJ200();
  else if (ecal_material == 18)
    EcalMaterial = MyMaterials::BGO();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << ecal_material << G4endl;
    exit(-1);
  }
  G4cout << "ECAL material: " << EcalMaterial << G4endl;

  /************************************************************************************/
  ScintiMaterial = NULL;
  if (scinti_material == 1)
    ScintiMaterial = MyMaterials::Quartz();
  else if (scinti_material == 2)
    ScintiMaterial = MyMaterials::SiO2();
  else if (scinti_material == 3)
    ScintiMaterial = MyMaterials::SiO2_Ce();
  else if (scinti_material == 4)
    ScintiMaterial = MyMaterials::LuAG_Ce();
  else if (scinti_material == 5)
    ScintiMaterial = MyMaterials::YAG_Ce();
  else if (scinti_material == 6)
    ScintiMaterial = MyMaterials::LSO();
  else if (scinti_material == 7)
    ScintiMaterial = MyMaterials::LYSO();
  else if (scinti_material == 8)
    ScintiMaterial = MyMaterials::LuAG_undoped();
  else if (scinti_material == 9)
    ScintiMaterial = MyMaterials::GAGG_Ce();
  else if (scinti_material == 10)
    ScintiMaterial = MyMaterials::LuAG_Pr();
  else if (scinti_material == 11)
    ScintiMaterial = MyMaterials::PbF2();
  else if (scinti_material == 12)
    ScintiMaterial = MyMaterials::PlasticBC408();
  else if (scinti_material == 13)
    ScintiMaterial = MyMaterials::PlasticBC418();
  else if (scinti_material == 14)
    ScintiMaterial = MyMaterials::PWO();
  else if (scinti_material == 15)
    ScintiMaterial = MyMaterials::Acrylic();
  else if (scinti_material == 16)
    ScintiMaterial = MyMaterials::copper();
  else if (scinti_material == 17)
    ScintiMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << scinti_material << G4endl;
    exit(-1);
  }

  CherencMaterial = NULL;
  if (Cherenc_material == 1)
    CherencMaterial = MyMaterials::Quartz();
  else if (Cherenc_material == 2)
    CherencMaterial = MyMaterials::SiO2();
  else if (Cherenc_material == 3)
    CherencMaterial = MyMaterials::SiO2_Ce();
  else if (Cherenc_material == 4)
    CherencMaterial = MyMaterials::LuAG_Ce();
  else if (Cherenc_material == 5)
    CherencMaterial = MyMaterials::YAG_Ce();
  else if (Cherenc_material == 6)
    CherencMaterial = MyMaterials::LSO();
  else if (Cherenc_material == 7)
    CherencMaterial = MyMaterials::LYSO();
  else if (Cherenc_material == 8)
    CherencMaterial = MyMaterials::LuAG_undoped();
  else if (Cherenc_material == 9)
    CherencMaterial = MyMaterials::GAGG_Ce();
  else if (Cherenc_material == 10)
    CherencMaterial = MyMaterials::LuAG_Pr();
  else if (Cherenc_material == 11)
    CherencMaterial = MyMaterials::PbF2();
  else if (Cherenc_material == 12)
    CherencMaterial = MyMaterials::PlasticBC408();
  else if (Cherenc_material == 13)
    CherencMaterial = MyMaterials::PlasticBC418();
  else if (Cherenc_material == 14)
    CherencMaterial = MyMaterials::PWO();
  else if (Cherenc_material == 15)
    CherencMaterial = MyMaterials::Acrylic();
  else if (Cherenc_material == 16)
    CherencMaterial = MyMaterials::copper();
  else if (Cherenc_material == 17)
    CherencMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Cherenc_material << G4endl;
    exit(-1);
  }

  CherenpMaterial = NULL;
  if (Cherenp_material == 1)
    CherenpMaterial = MyMaterials::Quartz();
  else if (Cherenp_material == 2)
    CherenpMaterial = MyMaterials::SiO2();
  else if (Cherenp_material == 3)
    CherenpMaterial = MyMaterials::SiO2_Ce();
  else if (Cherenp_material == 4)
    CherenpMaterial = MyMaterials::LuAG_Ce();
  else if (Cherenp_material == 5)
    CherenpMaterial = MyMaterials::YAG_Ce();
  else if (Cherenp_material == 6)
    CherenpMaterial = MyMaterials::LSO();
  else if (Cherenp_material == 7)
    CherenpMaterial = MyMaterials::LYSO();
  else if (Cherenp_material == 8)
    CherenpMaterial = MyMaterials::LuAG_undoped();
  else if (Cherenp_material == 9)
    CherenpMaterial = MyMaterials::GAGG_Ce();
  else if (Cherenp_material == 10)
    CherenpMaterial = MyMaterials::LuAG_Pr();
  else if (Cherenp_material == 11)
    CherenpMaterial = MyMaterials::PbF2();
  else if (Cherenp_material == 12)
    CherenpMaterial = MyMaterials::PlasticBC408();
  else if (Cherenp_material == 13)
    CherenpMaterial = MyMaterials::PlasticBC418();
  else if (Cherenp_material == 14)
    CherenpMaterial = MyMaterials::PWO();
  else if (Cherenp_material == 15)
    CherenpMaterial = MyMaterials::Acrylic();
  else if (Cherenp_material == 16)
    CherenpMaterial = MyMaterials::copper();
  else if (Cherenp_material == 17)
    CherenpMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Cherenp_material << G4endl;
    exit(-1);
  }

  WrapMaterial = NULL;
  if (wrap_material == 1)
    WrapMaterial = MyMaterials::Quartz();
  else if (wrap_material == 2)
    WrapMaterial = MyMaterials::SiO2();
  else if (wrap_material == 3)
    WrapMaterial = MyMaterials::SiO2_Ce();
  else if (wrap_material == 4)
    WrapMaterial = MyMaterials::LuAG_Ce();
  else if (wrap_material == 5)
    WrapMaterial = MyMaterials::YAG_Ce();
  else if (wrap_material == 6)
    WrapMaterial = MyMaterials::LSO();
  else if (wrap_material == 7)
    WrapMaterial = MyMaterials::LYSO();
  else if (wrap_material == 8)
    WrapMaterial = MyMaterials::LuAG_undoped();
  else if (wrap_material == 9)
    WrapMaterial = MyMaterials::GAGG_Ce();
  else if (wrap_material == 10)
    WrapMaterial = MyMaterials::LuAG_Pr();
  else if (wrap_material == 11)
    WrapMaterial = MyMaterials::PbF2();
  else if (wrap_material == 12)
    WrapMaterial = MyMaterials::PlasticBC408();
  else if (wrap_material == 13)
    WrapMaterial = MyMaterials::PlasticBC418();
  else if (wrap_material == 14)
    WrapMaterial = MyMaterials::PWO();
  else if (wrap_material == 15)
    WrapMaterial = MyMaterials::Acrylic();
  else if (wrap_material == 16)
    WrapMaterial = MyMaterials::copper();
  else if (wrap_material == 17)
    WrapMaterial = MyMaterials::Epoxy();
  else if (wrap_material == 18)
    WrapMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << wrap_material << G4endl;
    exit(-1);
  }

  /************************************************************************************/
  WindowMaterial = MyMaterials::silicone();
  PMTGapMaterial = MyMaterials::silicone();
  GaMaterial = NULL;
  if (gap_material == 1)
    GaMaterial = MyMaterials::Air();
  else if (gap_material == 2)
    GaMaterial = MyMaterials::OpticalGrease();
  else if (gap_material == 3)
    GaMaterial = MyMaterials::MeltMount168();
  else if (gap_material == 4)
    GaMaterial = MyMaterials::OpticalGrease155();
  else if (gap_material == 5)
    GaMaterial = MyMaterials::silicone();
  else if (gap_material == 6)
    GaMaterial = MyMaterials::PyrexGlass();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
    exit(-1);
  }
  G4cout << "Gap material: " << gap_material << G4endl;

  DeMaterial = NULL;
  if (det_material == 1)
    DeMaterial = MyMaterials::Silicon();
  else if (det_material == 2)
    DeMaterial = MyMaterials::Quartz();
  else if (det_material == 3)
    DeMaterial = MyMaterials::Air();
  else if (det_material == 4)
    DeMaterial = MyMaterials::Bialkali();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
    exit(-1);
  }
  G4cout << "Detector material: " << det_material << G4endl;

  //------------------
  // change properties

  if (core_absLength > 0)
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = {1. * eV, 10. * eV};
    G4double Absorption[nEntries_ABS] = {core_absLength * mm, core_absLength * mm};

    CoMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
    CoMaterial->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH", PhotonEnergy_ABS, Absorption, nEntries_ABS);
  }
  if (core_rIndex > 0)
  {
    const G4int nEntries_RI = 2;
    G4double PhotonEnergy_RI[nEntries_RI] = {1. * eV, 10. * eV};
    G4double RefractiveIndex[nEntries_RI] = {core_rIndex, core_rIndex};

    CoMaterial->GetMaterialPropertiesTable()->RemoveProperty("RINDEX");
    CoMaterial->GetMaterialPropertiesTable()->AddProperty("RINDEX", PhotonEnergy_RI, RefractiveIndex, nEntries_RI);
  }
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void DetectorConstruction::ConstructField()
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl;

  static G4TransportationManager *trMgr = G4TransportationManager::GetTransportationManager();

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager *globalFieldMgr = trMgr->GetFieldManager();

  if (!B_field_IsInitialized)
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector(0.0522 * B_field_intensity, 0.0522 * B_field_intensity, 0.9973 * B_field_intensity);

    B_field = new G4UniformMagField(fieldVector);
    globalFieldMgr->SetDetectorField(B_field);
    globalFieldMgr->CreateChordFinder(B_field);
    globalFieldMgr->GetChordFinder()->SetDeltaChord(0.005 * mm);
    B_field_IsInitialized = true;
  }

  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl;
  return;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit) && (maxStep > 0.))
    stepLimit->SetMaxAllowedStep(maxStep);
}

// //-----------------------------------------------
// //------------- Surface properties --------------
// //-----------------------------------------------

void DetectorConstruction::initializeSurface()
{
  fFiberWrapSurface = MakeS_wrap();
  fPMTSurface = MakeS_PMT();
  fPMTCaseSurface = MakeS_IdealWhiteSurface();

  static const unsigned nentries = 2;
  static double phoE[nentries] = {0.01 * eV, 10.0 * eV};
  double reflectivity[nentries] = {wrap_ref, wrap_ref};
  G4MaterialPropertiesTable *table = fFiberWrapSurface->GetMaterialPropertiesTable();
  if (table)
  {
    table->RemoveProperty("REFLECTIVITY");
    table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  }
  else
  {
    table = new G4MaterialPropertiesTable();
    table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
    fFiberWrapSurface->SetMaterialPropertiesTable(table);
  }
  fIdealPolishedOpSurface = MakeS_IdealPolished();
  /*
//perfect reflection
  static const unsigned nentriesc = 2;
  static double phoEc[nentries]   = {0.01*eV, 10.0*eV};
  double reflectivityc[nentries]  = {1,1};
  G4MaterialPropertiesTable* tablecrys = fIdealPolishedOpSurface->GetMaterialPropertiesTable();
  if( tablecrys ){
    tablecrys->RemoveProperty( "REFLECTIVITY" );
    tablecrys->AddProperty( "REFLECTIVITY", phoEc, reflectivityc, nentriesc );
  } else {
    tablecrys = new G4MaterialPropertiesTable();
    table->AddProperty( "REFLECTIVITY", phoEc, reflectivityc, nentriesc );
    fIdealPolishedOpSurface->SetMaterialPropertiesTable( tablecrys );
  }
*/

  fFilterSurface_ff = MakeS_IdealPolished();
  fFilterSurface_fr = MakeS_IdealPolished();

  const G4int nEntries_tran_u330 = 18;
  G4double PhotonEnergy_tran_u330[nEntries_tran_u330] = {4.103235582 * eV, 4.048912761 * eV, 3.944470896 * eV, 3.81331407 * eV, 3.750947295 * eV, 3.690598697 * eV, 3.603625797 * eV, 3.534216059 * eV, 3.456978466 * eV, 3.380554073 * eV, 3.309819592 * eV, 3.230572067 * eV, 3.185706353 * eV, 3.131340814 * eV, 3.087086539 * eV, 3.050146549 * eV, 2.992445212 * eV, 2.933127681 * eV};
  G4double transIndex_u330[nEntries_tran_u330] = {0.201372, 0.202705, 0.211043, 0.227125, 0.234102, 0.233987, 0.235942, 0.235798, 0.229958, 0.219856, 0.199831, 0.155664, 0.115833, 0.068881, 0.037554, 0.017592, 0.00466, 0.000935};

  const G4int nEntries_tran_ug5 = 25;
  G4double PhotonEnergy_tran_ug5[nEntries_tran_ug5] = {4.092143963 * eV, 4.034407045 * eV, 3.943981544 * eV, 3.825267479 * eV, 3.743879546 * eV, 3.62234533 * eV, 3.530100421 * eV, 3.414187953 * eV, 3.300875562 * eV, 3.233225815 * eV, 3.168293273 * eV, 3.137871012 * eV, 3.099604675 * eV, 3.074608111 * eV, 3.041899835 * eV, 3.001980276 * eV, 2.947821354 * eV, 2.873756177 * eV, 2.764363413 * eV, 2.697530944 * eV, 2.642988727 * eV, 2.564470256 * eV, 2.529030177 * eV, 2.498643446 * eV, 2.482374634 * eV};
  G4double transIndex_ug5[nEntries_tran_ug5] = {0.197435, 0.199462, 0.209962, 0.224666, 0.230184, 0.23489, 0.238226, 0.232947, 0.219825, 0.204682, 0.176017, 0.153858, 0.125979, 0.106681, 0.089491, 0.06658, 0.04575, 0.031245, 0.029397, 0.029832, 0.02319, 0.016426, 0.012699, 0.009703, 0.008199};

  G4MaterialPropertiesTable *Filtertable_u330 = new G4MaterialPropertiesTable();
  Filtertable_u330->AddProperty("REFLECTIVITY", PhotonEnergy_tran_u330, transIndex_u330, nEntries_tran_u330);
  G4MaterialPropertiesTable *Filtertable_ug5 = new G4MaterialPropertiesTable();
  Filtertable_ug5->AddProperty("REFLECTIVITY", PhotonEnergy_tran_ug5, transIndex_ug5, nEntries_tran_ug5);

  if (front_filter == -1)
  {
    std::cout << "No front Filter!" << std::endl;
  }
  else if (front_filter == 0)
  {
    fFilterSurface_ff->SetMaterialPropertiesTable(Filtertable_u330);
  }
  else if (front_filter == 1)
  {
    fFilterSurface_ff->SetMaterialPropertiesTable(Filtertable_ug5);
  }
  if (rear_filter == -1)
  {
    std::cout << "No rear Filter!" << std::endl;
  }
  else if (rear_filter == 0)
  {
    fFilterSurface_fr->SetMaterialPropertiesTable(Filtertable_u330);
  }
  else if (rear_filter == 1)
  {
    fFilterSurface_fr->SetMaterialPropertiesTable(Filtertable_ug5);
  }

  /*
  G4MaterialPropertiesTable* Filtertable_ff = fPMTSurface_ff->GetMaterialPropertiesTable();
  if( Filtertable ){
    Filtertable->RemoveProperty( "TRANSMITTANCE" );
    Filtertable->AddProperty("TRANSMITTANCE", PhotonEnergy_tran, transIndex, nEntries_tran);
  } else {
    Filtertable = new G4MaterialPropertiesTable();
    Filtertable->AddProperty("TRANSMITTANCE", PhotonEnergy_tran, transIndex, nEntries_tran);
    fFilterSurface->SetMaterialPropertiesTable( Filtertable );
  }
*/

  /*
  G4MaterialPropertiesTable* Filtertable = fFilterSurface->GetMaterialPropertiesTable();
  const G4int nEntries_tran = 18;
  G4double PhotonEnergy_tran[nEntries_tran] = {4.103235582*eV, 4.048912761*eV,3.944470896     *eV,3.81331407      *eV,3.750947295     *eV,3.690598697     *eV,3.603625797     *eV,3.534216059     *eV,3.456978466     *eV,3.380554073     *eV,3.309819592     *eV,3.230572067     *eV,3.185706353     *eV,3.131340814     *eV,3.087086539     *eV,3.050146549     *eV,2.992445212     *eV,2.933127681     *eV};

  G4double transIndex[nEntries_tran] = {0.201372        ,0.202705        ,0.211043        ,0.227125        ,0.234102        ,0.233987        ,0.235942        ,0.235798        ,0.229958        ,0.219856        ,0.199831        ,0.155664        ,0.115833        ,0.068881        ,0.037554        ,0.017592        ,0.00466 ,0.000935        };
  Filtertable->AddProperty("TRANSMITTANCE", PhotonEnergy_tran, transIndex, nEntries_tran);
  fFilterSurface->SetMaterialPropertiesTable( Filtertable );
*/
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

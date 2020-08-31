from pathlib import Path
import pickle
import unittest

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_almost_equal

from morfeus.morfeus import LocalForce, Sterimol, BuriedVolume, ConeAngle, SASA, Dispersion
from morfeus.io import read_xyz, read_gjf
from morfeus.helpers import convert_elements

CURR_DIR = Path(__file__).parent

class TestHelpers(unittest.TestCase):
    def test_convert_elements(self):
        """
        Test that elements can be converted back and forth from atomic numbers.
        """
        atomic_numbers = [1, 6, 7, 8, 9]
        atomic_symbols = ["H", "C", "N", "O", "F"]
        converted_numbers = convert_elements(atomic_numbers, output="symbols")
        converted_symbols = convert_elements(atomic_symbols)
        self.assertEqual(converted_numbers, atomic_symbols)
        self.assertEqual(converted_symbols, atomic_numbers)

class TestSterimol(unittest.TestCase):
    def setUp(self):
        self.sm_dir = CURR_DIR / "data/sterimol"

    def test_radii(self):
        radii_test = {
            "crc": {"L_value": 3.1999999999999997,
                    "L_value_uncorrected": 2.8,
                    "B_1_value": 1.7,
                    "B_5_value": 2.1358087170087674,
                    "L": np.array([2.8, 0., 0.]),
                    "B_1": np.array([0., 1.67667259, 0.28065817]),
                    "B_5": np.array([0., 0.55563877, 2.06226682]),
                    },    
            "bondi": {"L_value": 3.1999999999999997,
                    "L_value_uncorrected": 2.8,
                    "B_1_value": 1.7167757593031965,
                    "B_5_value": 2.2358087170087675,
                    "L": np.array([2.8, 0., 0.]),
                    "B_1": np.array([0., -1.21473827, 1.21314878]),
                    "B_5": np.array([0., 0.58165415, 2.15882354]),},
            }        
        for radii_type, results in radii_test.items():
            with self.subTest(radii_type=radii_type):            
                elements, coordinates = read_gjf(self.sm_dir / "me.gjf")
                sterimol = Sterimol(elements, coordinates, 1, 2, radii_type=radii_type)

                for attribute, value in results.items():
                    if isinstance(value, np.ndarray):
                        assert_array_almost_equal(getattr(sterimol, attribute), value)
                    else:
                        self.assertAlmostEqual(getattr(sterimol, attribute), value)

    def test_rot_vectors(self):
        rot_vectors_test = {
            3600: {"L_value": 3.1999999999999997,
                   "L_value_uncorrected": 2.8,
                   "B_1_value": 1.7,
                   "B_5_value": 2.1358087170087674,
                   "L": np.array([2.8, 0., 0.]),
                   "B_1": np.array([0., 1.67667259, 0.28065817]),
                   "B_5": np.array([0., 0.55563877, 2.06226682]),
                   },
            10: {"L_value": 3.1999999999999997,
                 "L_value_uncorrected": 2.8,
                 "B_1_value": 1.83,
                 "B_5_value": 2.132038059483839,
                 "L": np.array([2.8, 0., 0.]),
                 "B_1": np.array([0., 1.83, 0.]),
                 "B_5": np.array([0., -2.00346043, -0.72919996]),
                 },       
            }
        for n_rot_vectors, results in rot_vectors_test.items():
            with self.subTest(n_rot_vectors=n_rot_vectors):            
                elements, coordinates = read_gjf(self.sm_dir / "me.gjf")
                sterimol = Sterimol(elements, coordinates, 1, 2, n_rot_vectors=n_rot_vectors)

                for attribute, value in results.items():
                    if isinstance(value, np.ndarray):
                        assert_array_almost_equal(getattr(sterimol, attribute), value, decimal=5)
                    else:
                        self.assertAlmostEqual(getattr(sterimol, attribute), value, decimal=5)                        

class TestLocalForce(unittest.TestCase):
    def setUp(self):
        self.xtb_dir = CURR_DIR / "data/local_force/xtb"
        self.unimovib_dir = CURR_DIR / "data/local_force/unimovib"
        self.gaussian_dir = CURR_DIR / "data/local_force/gaussian"

    def test_gaussian_compliance(self):
        with open(self.gaussian_dir / "compliance/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.gaussian_dir / "compliance/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        lf = LocalForce()
        lf.load_file(self.gaussian_dir / "compliance/freq.fchk", "gaussian", "fchk")
        lf.compute_compliance()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)

    def test_gaussian_fchk(self):
        with open(self.gaussian_dir / "fchk/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.gaussian_dir / "fchk/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        lf = LocalForce()
        lf.load_file(self.gaussian_dir / "fchk/freq.fchk", "gaussian", "fchk")
        lf.normal_mode_analysis()
        lf.compute_local()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)

    def test_gaussian_local(self):
        with open(self.gaussian_dir / "local/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.gaussian_dir / "local/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        lf = LocalForce()
        lf.load_file(self.gaussian_dir / "local/freq.log", "gaussian", "log")
        lf.compute_local()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)

    def test_gaussian_hp(self):
        with open(self.gaussian_dir / "hp/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.gaussian_dir / "hp/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        lf = LocalForce()
        lf.load_file(self.gaussian_dir / "hp/freq.log", "gaussian", "log")
        lf.compute_local()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)    

    def test_misc(self):
        elements, coordinates = read_xyz(self.xtb_dir / "hessian/xtb.xyz")       
        lf = LocalForce(elements, coordinates)
        lf.load_file(self.xtb_dir / "hessian/hessian", "xtb", "hessian")
        lf.normal_mode_analysis()
        lf.add_internal_coordinate([1, 2])
        lf.add_internal_coordinate([2, 1, 3])
        lf.add_internal_coordinate([1, 2, 3, 4])
        lf.compute_local()
        lf.compute_frequencies()

        assert_almost_equal(lf.get_local_force_constant([1, 2]),
                            5.190222259808879, decimal=4)
        assert_almost_equal(lf.get_local_force_constant([2, 1, 3]),
                            2.3714906326401217, decimal=4)
        assert_almost_equal(lf.get_local_force_constant([1, 2, 3, 4]),
                            3.5142429141308917, decimal=4)
        assert_almost_equal(lf.get_local_frequency([1, 2]),
                            3078.130379468432, decimal=2)
        assert_almost_equal(lf.get_local_frequency([2, 1, 3]),
                            1457.3466611738797, decimal=2)
        assert_almost_equal(lf.get_local_frequency([1, 2, 3, 4]),
                            1457.9180060719998, decimal=2)

    def test_unimovib_local(self):
        with open(self.unimovib_dir / "local/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.unimovib_dir / "local/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        lf = LocalForce()
        lf.load_file(self.unimovib_dir / "local/localmode.dat", "unimovib", "local")
        lf.normal_mode_analysis()
        lf.detect_bonds()
        lf.compute_local()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)

    def test_unimovib_log(self):
        with open(self.unimovib_dir / "log/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.unimovib_dir / "log/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        lf = LocalForce()
        lf.load_file(self.unimovib_dir / "log/job.out", "unimovib", "log")
        lf.detect_bonds()
        lf.compute_local()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)

    def test_unimovib_umv(self):
        with open(self.unimovib_dir / "umv/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.unimovib_dir / "umv/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        lf = LocalForce()
        lf.load_file(self.unimovib_dir / "umv/job.umv", "unimovib", "umv")
        lf.normal_mode_analysis()
        lf.detect_bonds()
        lf.compute_local()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)

    def test_xtb_hessian(self):
        with open(self.xtb_dir / "hessian/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.xtb_dir / "hessian/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        elements, coordinates = read_xyz(self.xtb_dir / "hessian/xtb.xyz")       
        lf = LocalForce(elements, coordinates)
        lf.load_file(self.xtb_dir / "hessian/hessian", "xtb", "hessian")
        lf.normal_mode_analysis()
        lf.detect_bonds()
        lf.compute_local()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)           

    def test_xtb_compliance(self):
        elements, coordinates = read_xyz(self.xtb_dir / "compliance/xtb.xyz")       
        lf = LocalForce(elements, coordinates)
        lf.load_file(self.xtb_dir / "compliance/hessian", "xtb", "hessian")
        lf.normal_mode_analysis()
        lf.detect_bonds()
        lf.compute_compliance()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants,
                                  np.array([5.56450853, 21.16851423]),
                                  decimal=4)
        assert_array_almost_equal(lf.local_frequencies,
                                  np.array([3187.18622363, 2357.92562845]),
                                  decimal=2)

        lf.normal_mode_analysis(save_hessian=True)
        lf.compute_compliance()
        lf.compute_frequencies()
        assert_array_almost_equal(lf.local_force_constants,
                                  np.array([5.74915802, 21.56202739]),
                                  decimal=4)
        assert_array_almost_equal(lf.local_frequencies,
                                  np.array([3239.63554963, 2379.74109906]),
                                  decimal=2)        

    def test_xtb_normal_modes(self):
        with open(self.xtb_dir / "normal_modes/force_constants.pickle", "rb") as file:
            ref_force_constants = pickle.load(file)
        with open(self.xtb_dir / "normal_modes/frequencies.pickle", "rb") as file:
            ref_frequencies = pickle.load(file) 

        elements, coordinates = read_xyz(self.xtb_dir / "normal_modes/xtb.xyz")       
        lf = LocalForce(elements, coordinates)
        lf.load_file(self.xtb_dir / "normal_modes/xtb_normalmodes", "xtb", "normal_modes")
        lf.detect_bonds()
        lf.compute_local()
        lf.compute_frequencies()
        
        assert_array_almost_equal(lf.local_force_constants, ref_force_constants, decimal=4)
        assert_array_almost_equal(lf.local_frequencies, ref_frequencies, decimal=2)

class TestBuriedVolume(unittest.TestCase):
    def setUp(self):
        self.bv_dir = CURR_DIR / "data/buried_volume"
    
    def test_parameters(self):
        tests = {
            "central_atom": [
                (1, 0.2962110976518822),
                (2, 0.175666231830041),
                ],
            "excluded_atoms": [
                ([1, 2, 3, 4, 5, 6, 7], 0.2962110976518822),
                ([1, 2, 3, 4], 0.5436020778978755),
                ],
            "include_hs": [
                (False, 0.2962110976518822),
                (True, 0.33532309914275066),
                ],
            "radius": [
                (3.5, 0.2962110976518822),
                (3.0, 0.28836412395709177),
                ],
            "radii_type": [
                ("bondi", 0.33532309914275066),
                ("crc", 0.3189643123369363),
                ],
            "radii_scale": [
                (1.17, 0.2962110976518822),
                (1.0, 0.2281611069698099),
                ],
            "density": [
                (0.001, 0.2962110976518822),
                (0.01, 0.29601505711318793)
                ],
            }
        for parameter, test in tests.items():
            with self.subTest(parameter=parameter, test=test):            
                elements, coordinates = read_xyz(self.bv_dir / "1.xyz")
                parameters = {
                    "central_atom": 1,
                    "excluded_atoms": [1, 2, 3, 4, 5, 6, 7],
                    "include_hs": False,
                    "radius": 3.5,
                    "radii_type": "bondi",
                    "radii_scale": 1.17,
                    "density": 0.001,
                    }
                for value, result in test:
                    if parameter == "radii_type":
                        parameters["include_hs"] = True
                    parameters[parameter] = value
                    bv = BuriedVolume(
                        elements,
                        coordinates,
                        parameters["central_atom"],
                        exclude_list=parameters["exclude_list"],
                        include_hs=parameters["include_hs"],
                        radius=parameters["radius"],
                        radii_type=parameters["radii_type"],
                        radii_scale=parameters["radii_scale"],
                        density=parameters["density"],
                        )

                    self.assertAlmostEqual(bv.buried_volume, result)

class TestConeAngle(unittest.TestCase):
    def setUp(self):
        self.ca_dir = CURR_DIR / "data/cone_angle"

    def test_radii_type(self):
        tests = [
            ("crc", 117.11012922937584, [6, 10, 13]),
            ("bondi", 120.4252396307234, [6, 10, 13]),
            ]

        for radii_type, cone_angle, tangent_atoms in tests:
            with self.subTest(radii_type=radii_type):            
                elements, coordinates = read_xyz(self.ca_dir / "PdPMe3.xyz")
                ca = ConeAngle(elements, coordinates, 1, radii_type=radii_type)
                self.assertAlmostEqual(ca.cone_angle, cone_angle)
                self.assertListEqual(ca.tangent_atoms, tangent_atoms)

#TODO Include atom areas in test
class TestSASA(unittest.TestCase):
    def setUp(self):
        self.sasa_dir = CURR_DIR / "data/sasa"
    
    def test_parameters(self):
        tests = {
            "radii_type": [
                ("crc", 288.3494102201136, 410.66140617653855),
                ("bondi", 293.0558760056855, 419.18215481025516),
                ],
            "probe_radius": [
                (1.4, 288.3494102201136, 410.66140617653855),
                (0.0, 147.06692524134047, 110.02054104931742),
                ],
            "density": [
                (0.01, 288.3494102201136, 410.66140617653855),
                (0.1, 287.78926555750263, 409.7685091969908),
                (0.001, 288.24559926866425, 410.50413130181255)
                ],
            }
        for parameter, test in tests.items():
            with self.subTest(parameter=parameter, test=test):            
                elements, coordinates = read_xyz(self.sasa_dir / "PdPMe3.xyz")
                parameters = {
                    "radii_type": "crc",
                    "probe_radius": 1.4,
                    "density": 0.01
                    }
                for value, area, volume in test:
                    parameters[parameter] = value
                    sasa = SASA(
                        elements,
                        coordinates,
                        probe_radius=parameters["probe_radius"],
                        radii_type=parameters["radii_type"],
                        density=parameters["density"],
                        )

                    self.assertAlmostEqual(sasa.area, area)    
                    self.assertAlmostEqual(sasa.volume, volume)


class TestDispersion(unittest.TestCase):
    def setUp(self):
        self.disp_dir = CURR_DIR / "data/dispersion"
    
    def test_d3(self):
        elements, coordinates = read_xyz(self.disp_dir / "acetonitrile.xyz")
        disp = Dispersion(elements, coordinates, compute_coefficients=False)
        disp.get_coefficients(self.disp_dir / "d3.out", model="d3")
        disp.calculate_p_int()

        self.assertAlmostEqual(disp.area, 88.46506645705581)
        self.assertAlmostEqual(disp.volume, 66.64120697173627)
        self.assertAlmostEqual(disp.p_int, 12.352745771163207)
        self.assertAlmostEqual(disp.p_max, 19.30185627069566)
        self.assertAlmostEqual(disp.p_min, 7.690653637045179)
    
    def test_d4(self):
        elements, coordinates = read_xyz(self.disp_dir / "acetonitrile.xyz")
        disp = Dispersion(elements, coordinates, compute_coefficients=False)
        disp.get_coefficients(self.disp_dir / "d4.out", model="d4")
        disp.calculate_p_int()

        self.assertAlmostEqual(disp.area, 88.46506645705581)
        self.assertAlmostEqual(disp.volume, 66.64120697173627)
        self.assertAlmostEqual(disp.p_int, 12.223352870737187)
        self.assertAlmostEqual(disp.p_max, 19.15625459877261)
        self.assertAlmostEqual(disp.p_min, 6.737379134171022)

    def test_cube(self):
        tests = {
            0.001: (84.9751049643602,
                    66.2841940738628,
                    12.243521875634714,
                    17.096730096567626,
                    7.302208606498548,
                    ),
            0.004: (63.1633493713322,
                    40.421168225017944,
                    23.943164950830617,
                    32.96594672918482,
                    14.277421613595518,
                    ),
            }
        for isodensity, results in tests.items():
            with self.subTest(isodensity=isodensity):
                area, volume, p_int, p_max, p_min = results
                elements, coordinates = read_xyz(self.disp_dir / "acetonitrile.xyz")
                disp = Dispersion(elements, coordinates, point_surface=False, compute_coefficients=False)
                disp.surface_from_cube(self.disp_dir / "density.cub", isodensity=isodensity)
                disp.get_coefficients()
                disp.calculate_p_int()
                
                self.assertAlmostEqual(disp.area, area)
                self.assertAlmostEqual(disp.volume, volume)
                self.assertAlmostEqual(disp.p_int, p_int)
                self.assertAlmostEqual(disp.p_max, p_max)
                self.assertAlmostEqual(disp.p_min, p_min)

    def test_multiwfn(self):
        elements, coordinates = read_xyz(self.disp_dir / "acetonitrile.xyz")
        disp = Dispersion(elements, coordinates, point_surface=False, compute_coefficients=False)
        disp.surface_from_multiwfn(self.disp_dir / "vtx.pdb")
        disp.get_coefficients()
        disp.calculate_p_int()
        
        self.assertAlmostEqual(disp.area, 84.76678423845655)
        self.assertAlmostEqual(disp.volume, 66.0597847021667)
        self.assertAlmostEqual(disp.p_int, 12.341064066451636)
        self.assertAlmostEqual(disp.p_max, 17.12992900588955)
        self.assertAlmostEqual(disp.p_min, 7.353049164413401)

    def test_id3(self):
        tests = {
            "radii_type": [
                ("rahm", 
                 88.46506645705581,
                 66.64120697173627,
                 12.349609564342398,
                 19.296840497225567,
                 7.686268240863259,
                 ),
                ("crc",
                 65.514045714077,
                 42.95413575935195,
                 22.036616400404434,
                 31.95023875363241,
                 15.071809649747706,
                 ),
                ],
            "density": [
                (0.1, 
                 88.46506645705581,
                 66.64120697173627,
                 12.349609564342398,
                 19.296840497225567,
                 7.686268240863259,
                 ),
                (0.01,
                 88.50798378916743,
                 66.59109504052978,
                 12.357057198246634,
                 20.103976085244426,
                 7.647239714738982,
                 ),
                ],
            #TODO add this test once the bug with excluded atoms is fixed
            #"excluded_atoms": [
            #    ([]),
            #    ([1, 2]),
            #    ],
            }
        for parameter, test in tests.items():
            with self.subTest(parameter=parameter, test=test):            
                elements, coordinates = read_xyz(self.disp_dir / "acetonitrile.xyz")
                parameters = {
                    "radii_type": "rahm",
                    "density": 0.1,
                    "excluded_atoms": [],
                    }
                for value, area, volume, p_int, p_max, p_min in test:
                    parameters[parameter] = value
                    disp = Dispersion(
                        elements,
                        coordinates,
                        radii_type=parameters["radii_type"],
                        density=parameters["density"],
                        excluded_atoms=parameters["excluded_atoms"],
                        )

                    self.assertAlmostEqual(disp.area, area)    
                    self.assertAlmostEqual(disp.volume, volume)     
                    self.assertAlmostEqual(disp.p_int, p_int)
                    self.assertAlmostEqual(disp.p_max, p_max)
                    self.assertAlmostEqual(disp.p_min, p_min)

if __name__ == '__main__':
    unittest.main()
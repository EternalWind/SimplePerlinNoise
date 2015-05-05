//#define PROFILING

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Runtime.CompilerServices;

using MathNet.Numerics.LinearAlgebra;

namespace SimplePerlinNoise
{
    public class ComplexNoise2D
    {
        private IEnumerable<Noise2D> m_Noises = null;
        private double m_TotalAmplitude = 0.0;

        internal ComplexNoise2D(IEnumerable<Noise2D> noises)
        {
            m_Noises = noises;
            
            foreach (var n in noises)
            {
                m_TotalAmplitude += n.Amplitude;
            }
        }

        /// <summary>
        /// Evaluate the noise at a normalised point (x, y).
        /// </summary>
        /// <param name="x">The normalised x coordinate</param>
        /// <param name="y">The normalised y coordinate</param>
        /// <returns>The value at the evaluated point</returns>
        public double EvalAt(double x, double y)
        {
            double result = 0.0;
            foreach (var n in m_Noises)
            {
                result += n.EvalAt(x, y);
            }

            //if (result > 1.0)
            //    result = 1.0;
            //else if (result < 0.0)
            //    result = 0.0;

            //return result;
            return (result + m_TotalAmplitude) / (m_TotalAmplitude * 2);
        }
    }

    public class Noise2D
    {
        private Vector<double>[,] m_Gradients;
        private int m_OffsetX;
        private int m_OffsetY;
        private int m_X;
        private int m_Y;
        private double m_Amplitude = 1.0;
        public double Amplitude
        {
            get
            {
                return m_Amplitude;
            }
        }

        internal Noise2D(Vector<double>[,] normalised_gradients, int x, int y, int x_offset, int y_offset, double amplitude)
        {
            m_Gradients = normalised_gradients;
            m_X = x;
            m_Y = y;
            m_OffsetX = x_offset;
            m_OffsetY = y_offset;
            m_Amplitude = amplitude;
        }

        /// <summary>
        /// Evaluate the noise at a normalised point (x, y).
        /// </summary>
        /// <param name="x">The normalised x coordinate</param>
        /// <param name="y">The normalised y coordinate</param>
        /// <returns>The value at the evaluated point</returns>
        public double EvalAt(double x, double y)
        {
            x = x * m_X;
            y = y * m_Y;

            int i = (int)Math.Floor(x);
            int j = (int)Math.Floor(y);

            var g00 = m_Gradients[i, j];
            var g10 = m_Gradients[i + 1, j];
            var g01 = m_Gradients[i, j + 1];
            var g11 = m_Gradients[i + 1, j + 1];

            double u = x - i;
            double v = y - j;

            double n00 = g00.DotProduct(Vector<double>.Build.DenseOfArray(new double[] { u, v }));
            double n10 = g10.DotProduct(Vector<double>.Build.DenseOfArray(new double[] { u - 1, v }));
            double n01 = g01.DotProduct(Vector<double>.Build.DenseOfArray(new double[] { u, v - 1 }));
            double n11 = g11.DotProduct(Vector<double>.Build.DenseOfArray(new double[] { u - 1, v - 1 }));

            double nx0 = Interpolate(n00, n10, u, F);
            double nx1 = Interpolate(n01, n11, u, F);

            double nxy = Math.Round(Interpolate(nx0, nx1, v, F), 5, MidpointRounding.AwayFromZero);
            /*
            Console.WriteLine("----------------------------------");
            Console.WriteLine("x: " + x);
            Console.WriteLine("y: " + y);
            Console.WriteLine("i: " + i);
            Console.WriteLine("j: " + j);
            Console.WriteLine();
            Console.WriteLine("g00: " + g00.ToString());
            Console.WriteLine("g10: " + g10.ToString());
            Console.WriteLine("g01: " + g01.ToString());
            Console.WriteLine("g11: " + g11.ToString());
            Console.WriteLine();
            Console.WriteLine("u: " + u);
            Console.WriteLine("v: " + v);
            Console.WriteLine("---------------------------------------------");
            */
            return nxy * m_Amplitude;
        }

        private double Interpolate(double from, double to, double percent, Func<double, double> f)
        {
            double f_out = f(percent);

            return from * (1 - f_out) + to * f_out;
        }

        private double F(double t)
        {
            return (double)(6.0 * Math.Pow(t, 5.0) - 15.0 * Math.Pow(t, 4.0) + 10 * Math.Pow(t, 3.0));
        }
    }

    public class NoiseGenerator
    {
        private static int[] PERMUTATION = new int[]
        {
            151,160,137,91,90,15,
            131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
            190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
            88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
            77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
            102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
            135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
            5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
            223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
            129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
            251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
            49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
            138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
            151,160,137,91,90,15,
            131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
            190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
            88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
            77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
            102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
            135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
            5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
            223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
            129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
            251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
            49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
            138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180, 151
        };

        private IList<Vector<double>> m_Gradients = null;
        public IList<double[]> Gradients
        {
            get
            {
                return new List<double[]>(m_Gradients.Select(t => t.ToArray()));
            }
        }

        public int Dimensions
        {
            get;
            private set;
        }

        private int permutatedSeed;
        public int Seed
        {
            get;
            private set;
        }

        private bool m_IsInitialized = false;

        public long TimeForCalculatingGradients { get; private set; }

        public long TimeForCalculatingDotProducts { get; private set; }

        public long TimeForCalculatingInterpolations { get; private set; }

        private Stopwatch stopWatch = new Stopwatch();

        public void Initialize(int seed)
        {
            if (!m_IsInitialized)
            {
                Seed = seed;
                permutatedSeed = PERMUTATION[seed % PERMUTATION.Length];

                m_IsInitialized = true;
            }
        }

        public void Initialize(int dimensions, int base_gradient_count, int seed)
        {
            if (!m_IsInitialized)
            {
                Dimensions = dimensions;
                Seed = seed;

                if (base_gradient_count <= 0)
                    throw new Exception("The base gradient count cannot be less than 1!");

                m_Gradients = CreateBaseGradients(base_gradient_count, seed);

                m_IsInitialized = true;
            }
        }

        public void Initialize(IList<Vector<double>> base_gradients, int seed)
        {
            if (!m_IsInitialized)
            {
                if (base_gradients.Count <= 0)
                    throw new Exception("The base gradient count cannot be less than 1!");

                Dimensions = base_gradients[0].Count;
                Seed = seed;

                if (base_gradients.Any(t => t.Count != Dimensions))
                    throw new Exception("The dimension for each base gradient needs to be the same!");

                m_Gradients = base_gradients;

                m_IsInitialized = true;
            }
        }

        private double[][] gs = new double[][]
        {
            new double[]{1.0, -1.0},
            new double[]{1.0, 0.0},
            new double[]{1.0, 1.0},

            new double[]{0.0, -1.0},
            new double[]{0.0, 0.0},
            new double[]{0.0, 1.0},

            new double[]{-1.0, 0.0},
            new double[]{-1.0, 1.0},
            new double[]{-1.0, -1.0},

            new double[]{1.0, -1.0},
            new double[]{1.0, 0.0},
            new double[]{1.0, 1.0},

            new double[]{0.0, -1.0},
            new double[]{0.0, 0.0},
            new double[]{0.0, 1.0},

            new double[]{-1.0, 0.0}
        };

        /// <summary>
        /// Generate a 2D perlin noise.
        /// </summary>
        /// <param name="x">The width of the noise grid. The number of cells in a row.</param>
        /// <param name="y">The height of the noise grid. The number of cells in a column.</param>
        /// <param name="x_offset">The starting x coordiate of the noise gird. The number of cells skipped from the origin along the x-axis.</param>
        /// <param name="y_offset">The starting y coordinate of the noise grid. The number of cells skipped from the origin along the y-axis.</param>
        /// <param name="amplitude">The scaling of the noise. The noise value ranges between -amplitude and amplitude.</param>
        /// <returns>The generated noise.</returns>
        public Noise2D Generate2D(int x, int y, int x_offset = 0, int y_offset = 0, double amplitude = 1.0)
        {
            if (!m_IsInitialized)
                throw new Exception("The generator needs to be initialized first!");

            if (Dimensions != 2)
                throw new Exception("The expected number of dimensions for output is different than the one used during initialization!");

            // Normally, there are n + 1 control points for n cells. However, we need 1 more control point for being able to evaluate at the ending boundary.
            // E.q. EvalAt(1.0, 1.0).
            int control_points_per_row = x + 2;
            int control_points_per_column = y + 2;

            Vector<double>[,] gradients = new Vector<double>[control_points_per_row, control_points_per_column];
            for (var i = 0; i < control_points_per_row; ++i)
                for (var j = 0; j < control_points_per_column; ++j)
                {
                    int index = GetHashForCoordinate(i + x_offset, j + y_offset) % m_Gradients.Count;

                    gradients[i, j] = m_Gradients[index];

                    //Console.WriteLine("(" + i + ", " + j + "): " + gradients[i, j].ToString());
                }

            Noise2D noise = new Noise2D(gradients, x, y, x_offset, y_offset, amplitude);

            return noise;
        }

        /// <summary>
        /// Generate a 2D perlin noise with fractal brownian motion.
        /// </summary>
        /// <param name="iterations">The number of iterations to apply. The number of noises to add together.</param>
        /// <param name="persistence">The persistence of fractal brownian motion. The attenuation speed of the amplitude for each iteration.</param>
        /// <param name="base_x">The width of the base noise grid. The number of cells in a row.</param>
        /// <param name="base_y">The height of the base noise grid. The number of cells in a column.</param>
        /// <param name="base_x_offset">The starting x coordiate of the base noise gird. The number of cells skipped from the origin along the x-axis.</param>
        /// <param name="base_y_offset">The starting y coordinate of the base noise grid. The number of cells skipped from the origin along the y-axis.</param>
        /// <returns>The generated noise.</returns>
        public ComplexNoise2D GenerateComplex2D(int iterations, double persistence, int base_x, int base_y, int base_x_offset = 0, int base_y_offset = 0)
        {
            if (iterations < 1)
                throw new Exception("You need to at least iterate once!");

            Noise2D[] noises = new Noise2D[iterations];

            for (int i = 0; i < iterations; ++i)
            {
                int frequency = (int)(Math.Pow(2.0, i) + 0.5);
                double amplitude = (double)Math.Pow(persistence, i);

                int x = frequency * base_x;
                int y = frequency * base_y;

                int x_offset = (int)((decimal)base_x_offset / base_x * x + 0.5m);
                int y_offset = (int)((decimal)base_y_offset / base_y * y + 0.5m);

                noises[i] = Generate2D(x, y, x_offset, y_offset, amplitude);

                //Console.WriteLine("----------------------");
            }

            return new ComplexNoise2D(noises);
        }

        public double Generate2D(double x, double y)
        {
#if PROFILING
            stopWatch.Restart();
#endif

            int i = (int)x;
            int j = (int)y;

            int mod_i = i & 255;
            int mod_j = j & 255;

            var permutation_length = PERMUTATION.Length;
            var gs_length = gs.Length - 1;

            var idx_g00 = PERMUTATION[permutatedSeed + mod_i] + mod_j;
            var idx_g01 = PERMUTATION[idx_g00 + 1] & gs_length;
            idx_g00 = PERMUTATION[idx_g00] & gs_length;

            var idx_g10 = PERMUTATION[permutatedSeed + mod_i + 1] + mod_j;
            var idx_g11 = PERMUTATION[idx_g10 + 1] & gs_length;
            idx_g10 = PERMUTATION[idx_g10] & gs_length;

#if PROFILING
            TimeForCalculatingGradients += stopWatch.ElapsedTicks;
            stopWatch.Restart();
#endif

            double u = x - i;
            double v = y - j;

            double n00 = CalculateWeightedGradient(idx_g00, u, v);
            double n10 = CalculateWeightedGradient(idx_g10, u - 1, v);
            double n01 = CalculateWeightedGradient(idx_g01, u, v - 1);
            double n11 = CalculateWeightedGradient(idx_g11, u - 1, v - 1);

#if PROFILING
            TimeForCalculatingDotProducts += stopWatch.ElapsedTicks;
            stopWatch.Restart();
#endif

            double fx = F(u);
            double fy = F(v);

            double nx0 = n00 + fx * (n10 - n00);
            double nx1 = n01 + fx * (n11 - n01);

            double nxy = nx0 + fy * (nx1 - nx0);

#if PROFILING
            TimeForCalculatingInterpolations += stopWatch.ElapsedTicks;
#endif

            return nxy;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private double CalculateWeightedGradient(int hash, double u, double v)
        {
            switch (hash)
            {
                case 0:
                    return u - v;
                case 1:
                    return u;
                case 2:
                    return u + v;
                case 3:
                    return -v;
                case 4:
                    return 0.0;
                case 5:
                    return v;
                case 6:
                    return -u;
                case 7:
                    return v - u;
                case 8:
                    return -u - v;
                case 9:
                    return u - v;
                case 10:
                    return u;
                case 11:
                    return u + v;
                case 12:
                    return -v;
                case 13:
                    return 0.0;
                case 14:
                    return v;
                case 15:
                    return -u;
                default:
                    return 0.0;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private double Interpolate(double from, double to, double percent, Func<double, double> f)
        {
            double f_out = f(percent);

            return from * (1 - f_out) + to * f_out;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private double F(double t)
        {
            return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
        }

        private IList<Vector<double>> CreateBaseGradients(int base_gradient_count, int seed)
        {
            var r = new Random(seed);

            var base_gradients = new Vector<double>[base_gradient_count];
            for (int i = 0; i < base_gradient_count; ++i)
            {
                //do
                //{
                //    base_gradients[i] = gs[r.Next(0, 9)];
                //    //base_gradients[i] = Vector<double>.Build.Random(Dimensions, r.Next());
                //} while (base_gradients[i].L2Norm() == 0.0);

                base_gradients[i] = base_gradients[i].Normalize(2.0);
            }

            return base_gradients;
        }

        private int GetNaturalNumberFor(int integer)
        {
            if (integer < 0)
                return -2 * integer - 1;
            else
                return 2 * integer;
        }

        private int GetHashForCoordinate(int x, int y)
        {
            x = GetNaturalNumberFor(x);
            y = GetNaturalNumberFor(y);

            int max = Math.Max(x, y);
            if (x == max)
                return x * x + x + y;
            else
                return y * y + x;
        }
    }
}

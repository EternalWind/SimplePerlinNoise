#define DOUBLE_DEPTH

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Configuration;
using System.Drawing;
using System.Drawing.Imaging;
using System.Diagnostics;

using SimplePerlinNoise;

namespace Test
{
    class Program
    {
        static byte[] GetRawBytes(short[] data)
        {
            byte[] raw_bytes = new byte[data.Length * sizeof(short)];
            Buffer.BlockCopy(data, 0, raw_bytes, 0, raw_bytes.Length);

            return raw_bytes;
        }

        static byte[] GetRawBytes(byte[] data)
        {
            return data;
        }

        static double[] NormaliseNoise(double[] non_uniform_noise)
        {
            double max = non_uniform_noise[0];
            double min = non_uniform_noise[0];

            foreach (double n in non_uniform_noise)
            {
                if (n > max)
                    max = n;
                if (n < min)
                    min = n;
            }

            double range = max - min;

            for (int i = 0; i < non_uniform_noise.Length; ++i)
                non_uniform_noise[i] = (non_uniform_noise[i] - min) / range;

            return non_uniform_noise;
        }

        static double MiddlePass(double val, double min, double max, double min_displacing_val = 0.0, double max_displacing_val = 1.0)
        {
            if (val > max)
                return max_displacing_val;
            else if (val < min)
                return min_displacing_val;
            else
                return val;
        }

        static double MiddleCut(double val, double min, double max, double mid_displacing_val)
        {
            if (val > max || val < min)
                return val;
            else
                return mid_displacing_val;
        }

        static double Contrast(double val, double min, double max, double min_displacing_val = 0.0, double max_displacing_val = 1.0)
        {
            double result = MiddlePass(val, min, max, min_displacing_val, max_displacing_val);

            result = (result - min_displacing_val) / (max_displacing_val - min_displacing_val);

            return result;
        }

        static double[,] Perturb(double[,] vals, Noise2D noiseX, Noise2D noiseY, int displace_range)
        {
            int x = vals.GetLength(0);
            int y = vals.GetLength(1);

            double[,] temp_vals = new double[x, y];
            for (int i = 0; i < x; ++i)
            {
                for (int j = 0; j < y; ++j)
                {
                    double n1 = noiseX.EvalAt((double)i / (x - 1), (double)j / (y - 1));
                    double n2 = noiseY.EvalAt((double)i / (x - 1), (double)j / (y - 1));

                    int u = i + (int)(n1 * displace_range);
                    int v = j + (int)(n2 * displace_range);

                    if (u < 0)
                        u = 0;
                    else if (u >= x)
                        u = x - 1;

                    if (v < 0)
                        v = 0;
                    else if (v >= y)
                        v = y - 1;

                    temp_vals[i, j] = vals[u, v];
                }
            }

            return temp_vals;
        }

        static double[,] Erode(double[,] vals, double smoothness)
        {
            int x = vals.GetLength(0);
            int y = vals.GetLength(1);

            for (int i = 1; i < x - 1; ++i)
            {
                for (int j = 1; j < y - 1; ++j)
                {
                    double max = 0.0;
                    int max_u = 0;
                    int max_v = 0;

                    if (vals[i, j] < 0.5)
                    {
                        for (int u = -1; u <= 1; ++u)
                        {
                            for (int v = -1; v <= 1; ++v)
                            {
                                if (u != 0 || v != 0)
                                {
                                    var diff = vals[i, j] - vals[i + u, j + v];
                                    if (Math.Abs(diff) > max)
                                    {
                                        max = diff;
                                        max_u = u;
                                        max_v = v;
                                    }
                                }
                            }
                        }
                    }

                    if (max != 0 && Math.Abs(max) <= (smoothness / x))
                    {
                        double fix = max / 2;
                        vals[i, j] -= fix;
                        vals[i + max_u, j + max_v] += fix;
                    }
                }
            }

            return vals;
        }

        static double[,] Smoothen(double[,] vals)
        {
            int x = vals.GetLength(0);
            int y = vals.GetLength(1);

            for (int i = 1; i < x - 1; ++i)
            {
                for (int j = 1; j < y - 1; ++j)
                {
                    if (vals[i, j] < 0.6)
                    {
                        double total = 0.0;

                        for (int u = -1; u <= 1; ++u)
                        {
                            for (int v = -1; v <= 1; ++v)
                            {
                                total += vals[i + u, j + v];
                            }
                        }

                        vals[i, j] = total / 9;
                    }
                }
            }

            return vals;
        }

        static void Main(string[] args)
        {
            var stop_watch = new Stopwatch();

            stop_watch.Start();
            var settings = ConfigurationManager.AppSettings;

            int dimensions = 2;
            if (!Int32.TryParse(settings["Dimensions"], out dimensions))
                Console.WriteLine("Failed to read 'Dimensions' from the config file.");

            int base_gradient_count = 8;
            if (!Int32.TryParse(settings["Base_Gradient_Count"], out base_gradient_count))
                Console.WriteLine("Failed to read 'Base_Gradient_Count' from the config file.");

            int seed = 1;
            if (!Int32.TryParse(settings["Random_Seed"], out seed))
                Console.WriteLine("Failed to read 'Random_Seed' from the config file.");

            int iterations = 1;
            if (!Int32.TryParse(settings["Iterations"], out iterations))
                Console.WriteLine("Failed to read 'Iterations' from the config file.");

            double persistence = 108;
            if (!Double.TryParse(settings["Persistence"], out persistence))
                Console.WriteLine("Failed to read 'Persistence' from the config file.");

            int source_x = 256;
            if (!Int32.TryParse(settings["Source_X"], out source_x))
                Console.WriteLine("Failed to read 'Source_X' from the config file.");

            int source_y = 108;
            if (!Int32.TryParse(settings["Source_Y"], out source_y))
                Console.WriteLine("Failed to read 'Source_Y' from the config file.");

            int source_offset_x = 0;
            if (!Int32.TryParse(settings["Source_X_Offset"], out source_offset_x))
                Console.WriteLine("Failed to read 'Source_X_Offset' from the config file.");

            int source_offset_y = 0;
            if (!Int32.TryParse(settings["Source_Y_Offset"], out source_offset_y))
                Console.WriteLine("Failed to read 'Source_Y_Offset' from the config file.");

            int output_x = 2560;
            if (!Int32.TryParse(settings["Output_Resolution_X"], out output_x))
                Console.WriteLine("Failed to read 'Output_Resolution_X' from the config file.");

            int output_y = 1080;
            if (!Int32.TryParse(settings["Output_Resolution_Y"], out output_y))
                Console.WriteLine("Failed to read 'Output_Resolution_Y' from the config file.");

            string output_file = settings["Output_File_Name"];

            Console.WriteLine("Time for reading settings: " + stop_watch.ElapsedMilliseconds);
            stop_watch.Restart();

            var generator = new NoiseGenerator();
            generator.Initialize(seed);

            //var noise = generator.GenerateComplex2D(iterations, persistence, source_x, source_y, source_offset_x, source_offset_y);
            Console.WriteLine("Time for initialization: " + stop_watch.ElapsedMilliseconds);

            //var noise2 = generator.GenerateComplex2D(iterations, persistence, source_x, source_y, source_offset_x + 1, source_offset_y);

            //int offset_x = (int)((output_x - 1) * (double)1 / (source_x - 1));

            //for (int j = 0; j < output_y; ++j)
            //{
            //    double sampling_pos_y = (double)j / (output_y - 1);

            //    for (int i = 0; i < output_x - offset_x; ++i)
            //    {
            //        double sampling_pos_x2 = (double)(i) / (output_x - 1);
            //        double sampling_pos_x1 = sampling_pos_x2 + 1.0 / (source_x);

            //        var r1 = noise.EvalAt(sampling_pos_x1, sampling_pos_y);
            //        var r2 = noise2.EvalAt(sampling_pos_x2, sampling_pos_y);

            //        if (r1 != r2)
            //        {
            //            Console.WriteLine("Diff at (" + i + ", " + j + ")! r1 = " + r1 + " r2 = " + r2);
            //        }
            //    }
            //}

            //Bitmap output_preview = new Bitmap(output_x, output_y);
#if DOUBLE_DEPTH
            PixelFormat format = PixelFormat.Format16bppRgb565;
#else
            PixelFormat format = PixelFormat.Format8bppIndexed;
#endif

            Bitmap preview_output = new Bitmap(output_x, output_y, PixelFormat.Format32bppRgb);

            Bitmap output = new Bitmap(output_x, output_y, format);
            Rectangle rect = new Rectangle(0, 0, output_x, output_y);
            BitmapData bmp_data = output.LockBits(rect, ImageLockMode.ReadWrite, output.PixelFormat);
            IntPtr ptr = bmp_data.Scan0;

#if DOUBLE_DEPTH
            int bytes = Math.Abs(bmp_data.Stride) * output.Height / 2;
            short[] raw_memory = new short[bytes];
            int max_depth = 65535;
#else
            int bytes = Math.Abs(bmp_data.Stride) * output.Height;
            byte[] raw_memory = new byte[bytes];
            int max_depth = 255;
#endif

            System.Runtime.InteropServices.Marshal.Copy(ptr, raw_memory, 0, bytes);

            stop_watch.Restart();

            double[,] non_uniform_noise = new double[output_x, output_y];

            for (int j = 0; j < output_y; ++j)
            {
                double sampling_pos_y = (double)j / (output_y - 1);

                for (int i = 0; i < output_x; ++i)
                {
                    double sampling_pos_x = (double)i / (output_x - 1);
                    non_uniform_noise[i, j] = (generator.Generate2D((double)sampling_pos_x * 4, (double)sampling_pos_y * 4) + 1.0) / 2.0;
                    //non_uniform_noise[i, j] = Contrast(non_uniform_noise[i, j], 0.5, 0.7, 0.5, 0.7);
                    //non_uniform_noise[pixel_idx] = MiddleCut(non_uniform_noise[pixel_idx], 0.4, 0.7, 0.55);
                    //non_uniform_noise[pixel_idx] = MiddleCut(non_uniform_noise[pixel_idx], 0.2, 0.3, 0.25);
                }
            }

            Console.WriteLine("Time for generating noise: " + stop_watch.ElapsedMilliseconds + " (" + Stopwatch.Frequency + ")");
            Console.WriteLine("Generation profiling: ");
            Console.WriteLine("Time for calculating gradients ---- " + generator.TimeForCalculatingGradients);
            Console.WriteLine("Time for calculating dot products ---- " + generator.TimeForCalculatingDotProducts);
            Console.WriteLine("Time for calculating interpolations ---- " + generator.TimeForCalculatingInterpolations);

            stop_watch.Restart();

            //int frequency = 32;
            //non_uniform_noise = Perturb(non_uniform_noise, generator.Generate2D(frequency, frequency), generator.Generate2D(frequency, frequency, frequency, frequency), 32);

            //for (int i = 0; i < 10; ++i)
            //{
            //    non_uniform_noise = Erode(non_uniform_noise, 16.0);
            //}

            //non_uniform_noise = Smoothen(non_uniform_noise);

            Console.WriteLine("Time for post processing: " + stop_watch.ElapsedMilliseconds);
            stop_watch.Stop();

            for (int j = 0; j < output_y; ++j)
            {
                double sampling_pos_y = (double)j / output_y;

                for (int i = 0; i < output_x; ++i)
                {
                    double sampling_pos_x = (double)i / output_x;
                    int pixel_idx = j * output_y + i;

                    double normalised_sampling_result = non_uniform_noise[i, j];
                    int sampling_result = (int)(max_depth * normalised_sampling_result);
                    int preview_sampling_result = (int)(255 * normalised_sampling_result);

                    if (sampling_result > max_depth)
                        throw new Exception("The sampling result is larger than 1!");

                    preview_output.SetPixel(i, j, Color.FromArgb(preview_sampling_result, preview_sampling_result, preview_sampling_result));
#if DOUBLE_DEPTH
                    raw_memory[pixel_idx] = (short)sampling_result;
#else
                    raw_memory[pixel_idx] = (byte)sampling_result;
#endif
                }
            }

            System.Runtime.InteropServices.Marshal.Copy(raw_memory, 0, ptr, bytes);
            output.UnlockBits(bmp_data);

            output.Save(output_file + ".png", ImageFormat.Png);
            preview_output.Save(output_file + "_preview.png", ImageFormat.Png);

            using (System.IO.FileStream fs = new System.IO.FileStream(output_file + ".raw", System.IO.FileMode.Create))
            {
                byte[] raw_bytes = GetRawBytes(raw_memory);
                fs.Write(raw_bytes, 0, raw_bytes.Length);
            }

            Console.ReadKey();
        }
    }
}

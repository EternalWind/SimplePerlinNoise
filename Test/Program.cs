//#define DOUBLE_DEPTH

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Configuration;
using System.Drawing;
using System.Drawing.Imaging;

using SimplePerlinNoise;

namespace SimpleNoiseTest
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

        static double HighPass(double value, double threshold)
        {
            if (value < threshold)
                return 0.0f;
            else
                return value;
        }

        static void Main(string[] args)
        {
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

            double plane_threshold = 0.0;
            if (!Double.TryParse(settings["Plane_Threshold"], out plane_threshold))
                Console.WriteLine("Failed to read 'Plane_Threshold' from the config file.");

            var begin = DateTime.Now;

            var generator = new NoiseGenerator();
            generator.Initialize(dimensions, base_gradient_count, seed);

            var noise = generator.Generate2D(source_x, source_y, source_offset_x, source_offset_y);

            Bitmap preview_output = new Bitmap(output_x, output_y, PixelFormat.Format32bppRgb);
#if DOUBLE_DEPTH
            PixelFormat format = PixelFormat.Format16bppRgb565;
#else
            PixelFormat format = PixelFormat.Format8bppIndexed;
#endif

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
          
            for (int j = 0; j < output_y; ++j)
            {
                double sampling_pos_y = (double)j / output_y;

                for (int i = 0; i < output_x; ++i)
                {
                    double sampling_pos_x = (double)i / output_x;

                    double normalised_sampling_result = (noise.EvalAt(sampling_pos_x, sampling_pos_y) + 1) / 2;
                    normalised_sampling_result = HighPass(normalised_sampling_result, plane_threshold);

                    int sampling_result = (int)(max_depth * normalised_sampling_result);
                    int preview_sampling_result = (int)(255 * normalised_sampling_result);

                    if (sampling_result > max_depth)
                        throw new Exception("The sampling result is larger than 1!");

                    preview_output.SetPixel(i, j, Color.FromArgb(preview_sampling_result, preview_sampling_result, preview_sampling_result));

                    int pixel_idx = j * output_y + i;

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

            var end = DateTime.Now;

            Console.WriteLine("Done! Time used: " + (end - begin).TotalMilliseconds + "ms");
            Console.ReadKey();
        }
    }
}

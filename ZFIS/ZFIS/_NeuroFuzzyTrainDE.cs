namespace ZFIS
{
	class _NeuroFuzzyTrainDE : DE.Problem
	{
		_NeuroFuzzy NF;
		double[][] Inputs;
		double[][] Outputs;

		internal int parameter_number;
		//internal double[] parameters;
		internal double[] rmin;
		internal double[] rmax;
		internal double[,] LRanges;
		internal double[,] RRanges;

		bool AdaptiveRT = true; // Trainable Right Terms

		internal _NeuroFuzzyTrainDE(
			_NeuroFuzzy network,
			double[][] inputs, double[][] desoutputs,
			double[,] lranges, double[,] rranges
		)
		{
			NF = network;
			LRanges = lranges;
			RRanges = rranges;
			Inputs = inputs;
			Outputs = desoutputs;
			int npar;
			int p;
			int i;
			npar = 0; for (p = 0; p < NF.NLParams; p++) npar += NF.LParams[p].NTerms;
			parameter_number = 7 * npar; // inputs: 3 for A-part and 4 for B-part
			if (AdaptiveRT)
			{
				npar = 0; for (p = 0; p < NF.NRParams; p++) npar += NF.RParams[p].NTerms;
				parameter_number += 3 * npar; // outputs: 3 for A-part and 0 for B-part
			}
			//
			Trial = new double[parameter_number];
			NPars = parameter_number;
			//
			rmin = new double[parameter_number];
			rmax = new double[parameter_number];
			int k = 0;
			for (p = 0; p < NF.NLParams; p++) for (i = 0; i < NF.LParams[p].NTerms; i++)
				{
					rmin[k] = rmin[k + 1] = rmin[k + 2] = rmin[k + 3] = rmin[k + 4] = rmin[k + 5] = rmin[k + 6] = lranges[p, 0];
					rmax[k] = rmax[k + 1] = rmax[k + 2] = rmax[k + 3] = rmax[k + 4] = rmax[k + 5] = rmax[k + 6] = lranges[p, 1];
					k += 7;
				}
			if (AdaptiveRT)
			{
				for (p = 0; p < NF.NRParams; p++) for (i = 0; i < NF.RParams[p].NTerms; i++)
					{
						rmin[k] = rmin[k + 1] = rmin[k + 2] = rranges[p, 0];
						rmax[k] = rmax[k + 1] = rmax[k + 2] = rranges[p, 1];
						k += 3;
					}
			}
			k = 0;
			for (p = 0; p < NF.NLParams; p++) for (i = 0; i < NF.LParams[p].NTerms; i++)
				{
					Trial[k++] = NF.LParams[p].Terms[i].Value.AL;
					Trial[k++] = NF.LParams[p].Terms[i].Value.AM;
					Trial[k++] = NF.LParams[p].Terms[i].Value.AR;
					Trial[k++] = 0;
					Trial[k++] = NF.LParams[p].Terms[i].Value.BL;
					Trial[k++] = NF.LParams[p].Terms[i].Value.BR;
					Trial[k++] = 1;
				}
			if (AdaptiveRT)
			{
				for (p = 0; p < NF.NRParams; p++) for (i = 0; i < NF.RParams[p].NTerms; i++)
					{
						Trial[k++] = NF.RParams[p].Terms[i].Value.AL;
						Trial[k++] = NF.RParams[p].Terms[i].Value.AM;
						Trial[k++] = NF.RParams[p].Terms[i].Value.AR;
					}
			}
		}

		void LoadParameters(double[] parameters)
		{
			int p, i, k;
			// next line limits range
			//for(k=0; k<parameter_number; k++) if(parameters[k]<rmin[k]||parameters[k]>rmax[k]) return 100;
			k = 0;
			for (p = 0; p < NF.NLParams; p++) for (i = 0; i < NF.LParams[p].NTerms; i++)
				{
					double b1 = parameters[k + 3];
					double b2 = parameters[k + 4];
					double b3 = parameters[k + 5];
					double b4 = parameters[k + 6];
					b2 = b1 + System.Math.Abs(b2 - b1);
					b3 = b2 + System.Math.Abs(b3 - b2);
					b4 = b3 + System.Math.Abs(b4 - b3);
					double bl = 1;
					double br = 1;
					if (b4 > b1)
					{
						bl = (b2 - b1) / (b4 - b1);
						br = (b3 - b1) / (b4 - b1);
					}
					NF.LParams[p].Terms[i].Value =
						new _ZNumber3(
							parameters[k], parameters[k + 1], parameters[k + 2],
							bl, (bl + br) / 2, br
						);
					k += 7;
				}
			if (!AdaptiveRT) return;
			for (p = 0; p < NF.NRParams; p++) for (i = 0; i < NF.RParams[p].NTerms; i++)
				{
					NF.RParams[p].Terms[i].Value =
						new _ZNumber3(
							parameters[k], parameters[k + 1], parameters[k + 2],
							1, 1, 1
						);
					k += 3;
				}
		}

		override public double EvaluateCost(double[] parameters)
		{
			LoadParameters(parameters);
			return NF.MSE(Inputs, Outputs);
		}

		// V2
		//
		override public double EvaluateError(double[] parameters)
		{
			LoadParameters(parameters);
			//return NF.MaxError(Inputs, Outputs);
			//
			double error = 0;
			for (int p = 0; p < NF.NLParams; p++)
			{
				//
				error += System.Math.Max(0, LRanges[p, 0] - NF.LParams[p].Terms[0].Value.AM);
				error += System.Math.Max(0, NF.LParams[p].Terms[NF.LParams[p].NTerms - 1].Value.AM - LRanges[p, 1]);
				//
				error += System.Math.Max(0, NF.LParams[p].Terms[0].Value.AR - LRanges[p, 1]);
				error += System.Math.Max(0, LRanges[p, 0] - NF.LParams[p].Terms[NF.LParams[p].NTerms - 1].Value.AL);
				//
				for (int i = 0; i < NF.LParams[p].NTerms - 1; i++)
				{
					error += System.Math.Max(0, NF.LParams[p].Terms[i].Value.AM - NF.LParams[p].Terms[i + 1].Value.AM);
				}
			}
			if (!AdaptiveRT) return error;
			for (int p = 0; p < NF.NRParams; p++)
			{
				//
				error += System.Math.Max(0, RRanges[p, 0] - NF.RParams[p].Terms[0].Value.AM);
				error += System.Math.Max(0, NF.RParams[p].Terms[NF.RParams[p].NTerms - 1].Value.AM - RRanges[p, 1]);
				//
				error += System.Math.Max(0, NF.RParams[p].Terms[0].Value.AR - RRanges[p, 1]);
				error += System.Math.Max(0, RRanges[p, 0] - NF.RParams[p].Terms[NF.RParams[p].NTerms - 1].Value.AL);
				//
				for (int i = 0; i < NF.RParams[p].NTerms - 1; i++)
				{
					error += System.Math.Max(0, NF.RParams[p].Terms[i].Value.AM - NF.RParams[p].Terms[i + 1].Value.AM);
				}
				//
				for (int i = 0; i < NF.RParams[p].NTerms; i++)
				{
					error += System.Math.Max(
						0,
						System.Math.Min(
							NF.RParams[p].Terms[i].Value.AR - NF.RParams[p].Terms[i].Value.AM,
							NF.RParams[p].Terms[i].Value.AM - NF.RParams[p].Terms[i].Value.AL
						) -
						(RRanges[p, 1] - RRanges[p, 0])

					);
				}
			}
			return error;
			//
		}
		//

		internal double LoadBest()
		{
			LoadParameters(Solution);
			return NF.MaxError(Inputs, Outputs);
		}

		override public int Save()
		{
			return Save("denfopt.txt");
		}

		void WriteLine(System.IO.StreamWriter f, string str)
		{
			string[] ss = str.Split('\n');
			for (int i = 0; i < ss.Length; i++)
			{
				f.Write(ss[i]);
				f.Write(f.NewLine);
			}
		}

		void WriteLine(System.IO.StreamWriter f)
		{
			f.Write(f.NewLine);
		}

		public int Save(string resfilename)
		{
			System.IO.StreamWriter f = new System.IO.StreamWriter(resfilename, false, System.Text.Encoding.ASCII);
			LoadBest();
			//
			int p, i;
			p = 0;
			for (p = 0; p < NF.NLParams; p++)
			{
				WriteLine(f, NF.LParams[p].Label);
				//System.Console.Write(NF.LParams[p].Label+"\n");
				for (i = 0; i < NF.LParams[p].NTerms; i++)
				{
					WriteLine(
						f,
						//System.Console.WriteLine(
						NF.LParams[p].Terms[i].Label +
						" = " +
						NF.LParams[p].Terms[i].Value
					);
				}
				WriteLine(f);
				//System.Console.WriteLine();
			}
			for (p = 0; p < NF.NRParams; p++)
			{
				WriteLine(f, NF.RParams[p].Label);
				//System.Console.Write(NF.RParams[p].Label+"\n");
				for (i = 0; i < NF.RParams[p].NTerms; i++)
				{
					WriteLine(
						f,
						//System.Console.WriteLine(
						NF.RParams[p].Terms[i].Label +
						" = " +
						NF.RParams[p].Terms[i].Value
					);
				}
				WriteLine(f);
				//System.Console.WriteLine();
			}
			//
			WriteLine(f);
			WriteLine(f, "Max Error: " + NF.MaxError(Inputs, Outputs));
			WriteLine(f, "MSE: " + NF.MSE(Inputs, Outputs));
			WriteLine(f, "RMSE: " + NF.RMSE(Inputs, Outputs));
			f.Close();
			//
			return 1;
		}

	}
}

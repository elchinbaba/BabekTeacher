using System;

namespace Cluster
{
	public class Clustering
	{

		protected int NC_; // number of clusters
		protected int NX_; // number of data vectors
		protected int ND_; // dimension of data vector

		public int NC { get { return NC_; } }
		public int NX { get { return NX_; } }
		public int ND { get { return ND_; } }

		public int NClusters { get { return NC_; } }
		public int Dim { get { return ND_; } }

		protected double[][] C; // cluster centers C[cluster][dimension]
		protected double[][] X; // data vectors X[data point][dimension]

		public double[][] Centers { get { return C; } }
		public double[][] Data { get { return X; } }

		protected double[][] U; // data vectors memberships to clusters U[cluster][data point] - mind the difference with clusterz.cs

		public double[][] Membership { get { return U; } }

		protected static double M = 2; // real fuzzifier

		public double Fuzzifier { get { return M; } set { if (value > 1) M = value; } }

		protected double[] Xmin, Xmax;
		public double[] XMinBounds { get { return Xmin; } }
		public double[] XMaxBounds { get { return Xmax; } }

		public Clustering(double[,] x)
		{
			NX_ = x.GetLength(0);
			ND_ = x.GetLength(1);
			X = new double[NX_][];
			for (int i = 0; i < NX_; i++)
			{
				X[i] = new double[ND_];
				for (int k = 0; k < ND_; k++) X[i][k] = x[i, k];
			}
			ComputeXRanges();
			InitClusters(1);
		}

		void ComputeXRanges()
		{
			int i, k;
			Xmin = new double[ND_];
			Xmax = new double[ND_];
			for (k = 0; k < ND_; k++)
			{
				Xmin[k] = X[0][k];
				Xmax[k] = X[0][k];
				for (i = 0; i < NX_; i++)
				{
					if (Xmin[k] > X[i][k]) Xmin[k] = X[i][k];
					if (Xmax[k] < X[i][k]) Xmax[k] = X[i][k];
				}
			}
		}

		protected delegate double TDistanceF(double[] v1, double[] v2);

		protected TDistanceF SqDistanceF = SqEuclideanDistanceF;

		protected static double SqEuclideanDistanceF(double[] v1, double[] v2)
		{
			double dist = 0;
			int dim = v1.Length;
			for (int k = 0; k < dim; k++)
			{
				double diff = v1[k] - v2[k];
				dist += diff * diff;
			}
			return dist;
		}

		protected TDistanceF DistanceF = EuclideanDistanceF;

		protected static double EuclideanDistanceF(double[] v1, double[] v2)
		{
			return Math.Sqrt(SqEuclideanDistanceF(v1, v2));
		}

		public double[] CenterOfGravity(double[][] vectors, double[] weights)
		{
			return CenterOfGravity(new double[ND_], vectors, weights);
		}

		public double[] CenterOfGravity(double[][] vectors, double[] weights, double power)
		{
			return CenterOfGravity(new double[ND_], vectors, weights, power);
		}

		public double[] CenterOfGravity(double[] cog, double[][] vectors, double[] weights)
		{
			int i, k;
			for (k = 0; k < ND_; k++)
			{
				double s1 = 0;
				double s2 = 0;
				for (i = 0; i < NX_; i++)
				{
					s1 += vectors[i][k] * weights[i];
					s2 += weights[i];
				}
				cog[k] = s1 / s2;
			}
			return cog;
		}

		public double[] CenterOfGravity(double[] cog, double[][] vectors, double[] weights, double power)
		{
			int i, k;
			for (k = 0; k < ND_; k++)
			{
				double s1 = 0;
				double s2 = 0;
				for (i = 0; i < NX_; i++)
				{
					s1 += vectors[i][k] * Math.Pow(weights[i], power);
					s2 += Math.Pow(weights[i], power);
				}
				cog[k] = s1 / s2;
			}
			return cog;
		}

		virtual public void InitClusters()
		{
			InitClusters(NC_);
		}

		virtual public void InitClusters(int nv)
		{
			if (nv >= 1) NC_ = nv;
			Random rand = new Random();
			C = new double[NC_][];
			U = new double[NC_][];
			int j, k;
			for (j = 0; j < NC_; j++)
			{
				C[j] = new double[ND_];
				for (k = 0; k < ND_; k++) C[j][k] = Xmin[k] + (Xmax[k] - Xmin[k]) * rand.NextDouble();
				U[j] = new double[NX_];
			}
		}

		virtual public void InitClusters(double[][] v)
		{
			NC_ = v.Length;
			C = new double[NC_][];
			U = new double[NC_][];
			int j, k;
			for (j = 0; j < NC_; j++)
			{
				C[j] = new double[ND_];
				for (k = 0; k < ND_; k++) C[j][k] = v[j][k];
				U[j] = new double[NX_];
			}
		}

		public double GetMembership(int c, int i)
		{
			return U[c][i];
		}

		virtual public double GetMembership(int c, double[] x)
		{
			return 0;
		}

		virtual public double GetMembership(int c, int d, double x)
		{
			return 0;
		}

		public double MinDistanceBetweenClusterCenters()
		{
			double dist = DistanceF(C[0], C[1]);
			int c1, c2;
			for (c1 = 0; c1 < NC_; c1++)
				for (c2 = c1 + 1; c2 < NC_; c2++) dist = Math.Min(dist, DistanceF(C[c1], C[c2]));
			return dist;
		}

		public static bool VectorIsGreater(double[] x1, double[] x2, int n)
		{
			for (int k = 0; k < n; k++)
			{
				if (x1[k] > x2[k]) return true;
				if (x1[k] < x2[k]) return false;
			}
			return false; // equal
		}

		public static bool VectorIsGreater(double[] x1, double[] x2, int[] dimorder)
		{
			for (int k = 0; k < dimorder.Length; k++)
			{
				if (x1[dimorder[k]] > x2[dimorder[k]]) return true;
				if (x1[dimorder[k]] < x2[dimorder[k]]) return false;
			}
			return false; // equal
		}

		public static void ExchangeVectors(ref double[] x1, ref double[] x2)
		{
			double[] temp = x1;
			x1 = x2;
			x2 = temp;
		}

		public void SortClusters()
		{
			for (int gap = NC_ / 2; gap > 0; gap /= 2)
			{
				for (int i = gap; i < NC_; i++)
				{
					for (int j = i - gap; j >= 0 && VectorIsGreater(C[j], C[j + gap], ND_); j -= gap)
					{
						ExchangeVectors(ref C[j], ref C[j + gap]);
						ExchangeVectors(ref U[j], ref U[j + gap]);
					}
				}
			}
		}

		public void SortClusters(int[] dimorder)
		{
			for (int gap = NC_ / 2; gap > 0; gap /= 2)
			{
				for (int i = gap; i < NC_; i++)
				{
					for (int j = i - gap; j >= 0 && VectorIsGreater(C[j], C[j + gap], dimorder); j -= gap)
					{
						ExchangeVectors(ref C[j], ref C[j + gap]);
						ExchangeVectors(ref U[j], ref U[j + gap]);
					}
				}
			}
		}

		virtual public double GetValidity()
		{
			int i, j;
			double val1 = 0;
			for (i = 0; i < NX_; i++)
			{
				double umax = U[0][i];
				for (j = 1; j < NC_; j++) if (umax < U[j][i]) umax = U[j][i];
				val1 += umax;
			}
			double val2 = 0;
			for (j = 0; j < NC_ - 1; j++)
			{
				for (int j2 = j + 1; j2 < NC_; j2++)
					for (i = 0; i < NX_; i++) val2 += Math.Min(U[j][i], U[j2][i]);
			}
			val2 /= NC_ * (NC_ - 1) / 2;
			return (val1 - val2) / NX_;
		}

	}
}

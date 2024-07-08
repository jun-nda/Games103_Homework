using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class implicit_model : MonoBehaviour
{
	float 		t 		= 0.0333f;
	float 		mass	= 1;
	float		damping	= 0.99f;
	float 		rho		= 0.995f;
	float 		spring_k = 8000;
	int[] 		E;
	float[] 	L;
	Vector3[] 	V;

    // Start is called before the first frame update
    void Start()
    {
		Mesh mesh = GetComponent<MeshFilter> ().mesh;

		//Resize the mesh.
		int n=21;
		Vector3[] vertices  	= new Vector3[n*n];
		Vector2[] UV 	= new Vector2[n*n];
		int[] triangles	= new int[(n-1)*(n-1)*6];
		for(int j=0; j<n; j++)
		for(int i=0; i<n; i++)
		{
			vertices[j*n+i] = new Vector3(5-10.0f*i/(n-1), 0, 5-10.0f*j/(n-1)); //  简单算一下，就是在往里塞定点数据，看着有点玄学
			UV[j*n+i]=new Vector3(i/(n-1.0f), j/(n-1.0f));
		}
		int t=0;
		for(int j=0; j<n-1; j++)
		for(int i=0; i<n-1; i++)	
		{
			triangles[t*6+0]=j*n+i;
			triangles[t*6+1]=j*n+i+1;
			triangles[t*6+2]=(j+1)*n+i+1;
			triangles[t*6+3]=j*n+i;
			triangles[t*6+4]=(j+1)*n+i+1;
			triangles[t*6+5]=(j+1)*n+i;
			t++;
		} // 这样弄完是一个小方块面片
		mesh.vertices=vertices;
		mesh.triangles=triangles;
		mesh.uv = UV;
		mesh.RecalculateNormals();


		// Construct the original E
		int[] _E = new int[ triangles.Length * 2];
		for (int i = 0; i < triangles.Length; i += 3) 
		{
			// 一个E[k]里面只存了一个点，不像我们平时是存一个数对，
			// 他这个就是连着的两个索引的值就是一条边的两个点，说实话还真不太好理解呢
			_E[i * 2 + 0] = triangles[i + 0];
			_E[i * 2 + 1] = triangles[i + 1];
			_E[i * 2 + 2] = triangles[i + 1];
			_E[i * 2 + 3] = triangles[i + 2];
			_E[i * 2 + 4] = triangles[i + 2];
			_E[i * 2 + 5] = triangles[i + 0];
		}
		//Reorder the original edge list
		for (int i = 0; i < _E.Length; i += 2)
			if(_E[i] > _E[i + 1]) 
				Swap(ref _E[i], ref _E[i+1]);
		// 上面这个是把一条边的两个顶点的序号进行排序了，1，0这种就变换成0，1
		
		//Sort the original edge list using quicksort
		Quick_Sort (ref _E, 0, _E.Length / 2 - 1);

		// 初始化，说实话这代码写的真的很难让人看懂
		// 前面不是E[i]和E[i+1]是一条边嘛, 下面这个计算边的数量
		// 同时去重，如果E[i]和E[i+1] 和前面两个相同，则略过 
		// 这一步纯是在算数组长度呢...吐了
		int e_number = 0;
		for (int i = 0; i < _E.Length; i += 2)
			if (i == 0 || _E [i + 0] != _E [i - 2] || _E [i + 1] != _E [i - 1]) 
				e_number ++;
		
		E = new int[e_number * 2]; // 这里还是存顶点
		// 逻辑和上面一样，如果相邻的两个边的顶点序号都不一样，才塞到真正的边数组里
		for (int i = 0, e = 0; i < _E.Length; i += 2)
		{
			if (i == 0 || _E [i + 0] != _E [i - 2] || _E [i + 1] != _E [i - 1]) 
			{
				E[e * 2 + 0]=_E [i + 0];
				E[e * 2 + 1]=_E [i + 1];
				e ++;
			}	
		}
		
		// 初始化每条的长度
		L = new float[E.Length / 2];
		for (int e = 0; e < E.Length / 2; e++) 
		{
			int v0 = E[e * 2 + 0];
			int v1 = E[e * 2 + 1];
			L[e] = (vertices[v0] - vertices[v1]).magnitude;
		}

		// 初始化速度
		V = new Vector3[vertices.Length];
		for (int i = 0; i < V.Length; i++)
			V[i] = new Vector3 (0, 0, 0);
    }
    
    // 有点像二分法呢
    void Quick_Sort(ref int[] a, int l, int r)
	{
		int j;
		if(l < r)
		{
			j = Quick_Sort_Partition(ref a, l, r);
			Quick_Sort (ref a, l, j-1);
			Quick_Sort (ref a, j+1, r);
		}
	}

	int  Quick_Sort_Partition(ref int[] a, int l, int r)
	{
		int pivot_0, pivot_1, i, j;
		pivot_0 = a [l * 2 + 0];
		pivot_1 = a [l * 2 + 1];
		i = l;
		j = r + 1;
		while (true) 
		{
			do ++i; while( i<=r && (a[i*2]<pivot_0 || a[i*2]==pivot_0 && a[i*2+1]<=pivot_1));
			do --j; while(  a[j*2]>pivot_0 || a[j*2]==pivot_0 && a[j*2+1]> pivot_1);
			if(i>=j)	break;
			Swap(ref a[i*2], ref a[j*2]);
			Swap(ref a[i*2+1], ref a[j*2+1]);
		}
		Swap (ref a [l * 2 + 0], ref a [j * 2 + 0]);
		Swap (ref a [l * 2 + 1], ref a [j * 2 + 1]);
		return j;
	}

	void Swap(ref int a, ref int b)
	{
		int temp = a;
		a = b;
		b = temp;
	}

	void Collision_Handling()
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] vertices = mesh.vertices;
		
		//Handle colllision.

		mesh.vertices = vertices;
	}

	void Get_Gradient(Vector3[] X, Vector3[] X_hat, float t, Vector3[] G)
	{
		// 这里的个人理解，这个公式本身是隐式积分那个方程组算出来的，M是质量，重力有没有都行
		// 噢卧槽没有重力果然是不行的
		// Momentum and Gravity. 先算公式的前半部分
		for (int i = 0; i < X.Length; ++i)
		{
			G[i] = -mass * new Vector3(0, -9.8f, 0);
		}

		// Spring Force.
		for (int i = 0; i < E.Length; i += 2)
		{
			int xi = E[i];
			int xj = E[i + 1];
			G[xi] += spring_k * (1 - L[i / 2] / (X[xi] - X[xj]).magnitude) * (X[xi] - X[xj]);
			G[xj] -= spring_k * (1 - L[i / 2] / (X[xi] - X[xj]).magnitude) * (X[xi] - X[xj]);
		} // 为啥xi一定是正的，xj一定是负的嘞
	}

    // Update is called once per frame
	void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X 		= mesh.vertices;
		Vector3[] last_X 	= new Vector3[X.Length];
		Vector3[] X_hat 	= new Vector3[X.Length];
		Vector3[] G 		= new Vector3[X.Length];

		//Initial Setup.
		for(int i = 0; i < V.Length; ++i)
		{
			V[i] *= damping;
		}

		for (int i = 0; i < X_hat.Length; i++)
		{
			last_X[i] = X[i];
			// X[i]也要初始化一下，减少迭代次数？ 为什么 不初始化 牛顿的次数可能会多
			X_hat[i] = X[i] + t * V[i];
			//X[i] = (X_hat[i] + X[i]) / 2f;
		}
		//Debug.Log("帧数: " + Time.frameCount);
		for(int k = 0; k < 32; k++)
		{
			//Debug.Log("当前迭代次数: " + k);
			Get_Gradient(X, X_hat, t, G);
			
			//Update vertices by gradient.
			// 牛顿法，迭代32次
			for (int i = 0; i < X.Length; ++i)
			{
				if (i == 0 || i == 20) continue;
				
				// 求导数，然后求出一个x，然后再用这个x用作下一个循环的自变量
				Vector3 Xk = X[i] - G[i] / (1 / (t * t) + 4 * spring_k);
				X[i] = Xk;
			}
		}

		//Finishing.
		for (int i = 0; i < V.Length; i++)
		{
			V[i] += (X[i] - last_X[i]) / t;
		}
		
		mesh.vertices = X;

		Collision_Handling ();
		mesh.RecalculateNormals ();
	}
}

// 这里面的东西都是一些新奇玩意，但是基本没用 

// 梳排序
// 时间复杂度 O(n^2)
// 哪个闲的蛋疼的人发明的。。。 
template<class T>
void combsort11_32(int a_len, T *a) { 
	int g = a_len; 
	uint32_t t; 
	while (true) { 
		int flag = 1; 
		g = (((g * 10) / 13) > 1) ? ((g * 10) / 13) : 1; 
		if (g == 9 || g == 10)  g = 11;
		for (int i = 0; i + g < a_len; i++) { 
			if (a[i] > a[i + g]) { t = a[i], a[i] = a[i + g], a[i + g] = t; flag = 0; } 
		} 
		if (g == 1 && flag) break;
	} 
}
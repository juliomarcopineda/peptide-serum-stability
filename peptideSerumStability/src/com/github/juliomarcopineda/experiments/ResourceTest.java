package com.github.juliomarcopineda.experiments;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class ResourceTest {

	public static void main(String[] args) {
		InputStream test = ResourceTest.class.getResourceAsStream("/test.txt");
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(test));
		
		String line;
		try {
			while ((line = reader.readLine()) != null) {
				System.out.println(line);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

}

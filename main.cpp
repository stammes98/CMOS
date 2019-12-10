#include <stdio.h>
#include <string>
#include "CMOSUtils.h"

int main() {
  doLiveCamera(2, 500.0, "../CMOSsensor/CMOSsensor/cImages/GlassTarget/0DegTest3.FTS", true, "../CMOSsensor/CMOSsensor/cImages/GlassTarget/0DegBias", "../CMOSsensor/CMOSsensor/cImages/GlassTarget/0DegDark", "test.csv", true);
	printf("Camera complete\n");
  return 0;
}

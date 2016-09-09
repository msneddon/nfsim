
<?php


function countDownload ($email)
{
	mysql_connect('mysql.emonetnew.biology.yale.edu','nfwebreguser','nepo4nEFsim') or die('Error 1.');
	@mysql_select_db('nfsimwebregister') or die('error 2.');
	
	
	$getUserQuery="SELECT id,downLoadCount FROM nfRegUser WHERE email=\"$email\";";
	$answer=mysql_query($getUserQuery);
	$rows=mysql_numrows($answer);
	
	$i=0;
	while($i<$rows) {
		$id = mysql_result($answer,$i,"id");
		$count = mysql_result($answer,$i,"downLoadCount");
		$newCount = $count + 1;
		//echo("<br>count is at:$count, will be:$newCount<br>");
		
		// now we update the database
		$updateCommand = "UPDATE nfRegUser SET downloadCount=$newCount WHERE id=$id";
		if (!mysql_query($updateCommand))  {
  			die('internal database error 4.');
  		}

		
		$i++;	
	}
	
	
	
	
	//exit database
	mysql_close();
	
	//return
	return;
	
	
}

?>
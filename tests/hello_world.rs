// テストモジュールを定義します。
#[cfg(test)]
mod tests {
    // テスト対象の関数をインポートします。
    use super::*;

    // テスト関数を定義します。
    #[test]
    fn test_hello_world() {
        // テスト対象の関数を呼び出し、結果をアサートします。
        print!("{:?}",hello_world());
        assert_eq!(hello_world(), "Hello, world!");
    }
}

// テスト対象の関数を定義します。
fn hello_world() -> &'static str {
    "Hello, world!"
}
